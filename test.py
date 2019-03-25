import zlib
import gzip 
import time
import timeit
import os
import mmap
import io
import gc
import argparse
import psutil
import multiprocessing
import copy
import ctypes
import operator
import inspect
#from multiprocessing.process import current_process
#current_process()._config["tempdir"] = "/dev/shm"

print (time.strftime("%c"))

mypid=os.getpid()

py = psutil.Process(mypid)

parser = argparse.ArgumentParser(description='Dedupe.')
parser.add_argument('--threads','-t', type=int, default=int(multiprocessing.cpu_count()/2), help='threads to use')
parser.add_argument('filenames', nargs='+', help='files to process')

args = parser.parse_args()
cpus = args.threads
filenames = args.filenames


shared_arr = multiprocessing.RawArray(ctypes.c_int,[0]*(cpus+1)*4)
arr_index = 0

manager = multiprocessing.Manager()
data_manager = manager.dict()
error_manager = manager.dict()

lock_manager = manager.Lock()
locked = multiprocessing.RawValue(ctypes.c_bool, False)

def best_overlap(read1,qual1,read2,qual2):
	'''
	Finds the longest perfect overlap between two reads
	Read1
	---------------------------
				 ||||||||||||||
				 ---------------------------
									   Read2
	If no overlap is found it returns -1, otherwise it returns a tuple
	with the consensus of SEQ and QUAL
	'''
	#Defines the variable consensus defaulted to -1
	consensus=-1
	
	#Loop through al the possible subStrings of Read1 that could be the
	#overlap in Read 2
	#Offset always contains the position of the first byte of the query
	for r1_offset in range(len(read1)):
		#OBS!: r1_query = read1[r1_offset:]
		
		#Check if the subString of Read1 is the beginning of Read2
		if read2.startswith(read1[r1_offset:]):
			#Reset r1_offset to the first hit of the query in r1.
			#If there is one match only, it shouldn't change.
			r1_offset=read1.find(read1[r1_offset:])
			
			#Calculate r2_offset to the last hit of the query in r2.
			#If there is only one match its value should be 0.
			r2_offset=read2.rfind(read1[r1_offset:])
			
			#Check whether there are multiple hits for the query.
			#If so, abort due to ambiguity
			if (read1.count(read1[r1_offset:]) == 1) and (read2.count(read1[r1_offset:]) == 1):
				#As we have taken the leftmost hit in r1 and the rightmost hit in r2,
				#the sequence after the hit in r1 and before the hit in r2 will only
				#be the same in the case a true consensus is found.
				#OBS! Probably redundant
				if read1[r1_offset:] == read2[:r2_offset + len(read1)-r1_offset]:
					#print("Only one hit of the query, proceed")
					#OBS!: read1[r1_offset:] == len(read1)-r1_offset
					consensus = (read1+read2[len(read1)-r1_offset:], qual1+qual2[len(read1)-r1_offset:])
					#print(read1.ljust(len(consensus[0]),b"-"))
					#print(read2.rjust(len(consensus[0]),b"-"))
					break
			else:
				consensus=-1
	return consensus

def getChunk(sourceFile, i, n_chunks):
	'''
	Chunks sourceFile in n_chunks and returns chunk i.
	Only the chunk of interest is saved in memory.
	Returns a tuple (startByte, byteString) where startByte is the
	first byte of the chunk relative to the full-sized file and
	byteString contains a chunk of sourceFile from startByte.
	'''
	#Open the file as read-only and byte mode
	with open(sourceFile, "rb") as fh:
		#Obtain the size of the file
		size=os.fstat(fh.fileno()).st_size
		#Calculate the size of the chunks
		chunksize=round(size/n_chunks)
		#Locate the start postition for chunk i
		fh.seek(i*chunksize)
		
		assert i*chunksize == fh.tell(), "The pointer in the file and the start of the chunk relative to sourceFile differ. They should be equal"
		return (i*chunksize, fh.read(chunksize))

def safe_readline(byteStream):
	'''
	Returns the first line in the byteStream the it is not empty
	from the position that the pointer is in. It only returns an
	empty byteString when EOF is reached
	'''
	#Variable to store the line. Defaults to empty byteString
	line=b''
	
	#Read the byteStream until the first non-empty line is found or EOF is reached
	#while (line == b'' and byteStream.tell() < byteStream.getbuffer().nbytes):
	#while (line == b'' and byteStream.tell() < len(byteStream.getbuffer())):
	while (line == b'' and byteStream.tell() < len(byteStream.getvalue())):
		line=byteStream.readline().strip()
	
	return line



def dedupe(seq, qual, stored_read):
	#Check if reads are identical. Omit if so.
	if seq == stored_read[0]:
		#print("identical",)
		return (1,) #No need to "return (1,stored_read)" as the read is not going to be used
		#pass
	
	#Check if the new read is contained in the stored read. Omit if so.
	elif seq in stored_read[0]:
		#print("s_read is longer")
		return (2,stored_read)
		#pass
	
	#Check if the stored read is contained in the new read.
	elif stored_read[0] in seq:
		#print("seq is longer, storing")
		return (2,(seq,qual))
	
	#If reads were neither identical nor contained, try to find a perfect overlap
	else:
		#Try to find overlap:
		#read1
		#----------------
		#          ---------------
		#                    read2
		consensus = best_overlap(seq,qual,stored_read[0],stored_read[1])
		
		#Check if overlap was found
		if consensus == -1:
			#Try to find overlap:
			#read2
			#----------------
			#          ---------------
			#                    read1
			consensus=best_overlap(stored_read[0],stored_read[1],seq,qual)
		
		#Check if overlap was found
		if consensus == -1:
			return None
			
		#If overlap was found
		else:
			#Add consensus (contains SEQ and QUAL) to data dict
			#print("consensus found",)
			return (2,consensus)


def processReads(byteString, s_data, s_error, shared_arr_index = -1):
	private_counter = [0,0,0,0]

	#Counter contains an array("added", "identical", "extended", "error")
	if(shared_arr_index == -1):
		print("Counter not initialised. Exiting")
		exit()

	p_error={}
	p_data={}
	
	'''
	It adds the reads in byteString to the shared dictionary data
	if they are correct/extended/overlapped or to the error
	dictionary if they were not properly processed.
	Returns true if it finished succesfully or exits otherwise
	'''
	#Load the byteString into a byteStream to loop over it
	with io.BytesIO(byteString) as file:

		#Make sure we are at the beginning of the stream.
		file.seek(0)
		
		#Execute until we reach the EOF
		while file.tell() < len(byteString): 
			#Wait if a lock is acquired which means shared data and shared error are being updated
			while locked.value == True:
				time.sleep(0.1)
				#pass
			
			#TODO: Explanation of internal counter. Basically to avoid bottleneck by writting too often to a shared object
			if(sum(private_counter) % 1000 == 0):
				shared_arr[shared_arr_index] += private_counter[0]
				shared_arr[shared_arr_index + 1] += private_counter[1]
				shared_arr[shared_arr_index + 2] += private_counter[2]
				shared_arr[shared_arr_index + 3] += private_counter[3]
				private_counter = [0,0,0,0]

			#FASTQ code is organised in groups of 4 NON-EMPTY lines.
			#3rd line is useless and ignored.
			#safe_readline() makes sure we do not load empty lines on the way.
			name = safe_readline(file)
			seq = safe_readline(file)
			safe_readline(file)
			qual = safe_readline(file)
			
			#If for some reason the file has empty lines at the end of the file, safe_readline will return b''
			if (name==b'' or seq==b'' or qual==b''):
				assert name==seq==qual==b'', "Error: read not loaded correctly\nName: "+name.decode("UTF-8")+"\nRead: "+seq.decode("UTF-8")+"\n Qual: "+qual.decode("UTF-8")
				#print(os.getpid(),": EOF reached with newlines at the end of the file")
				break
				'''
				if name==seq==qual==b'':
					print(os.getpid(),": EOF reached with newlines at the end of the file")
				else:
					print("Error: read not loaded correctly")
					exit()
				'''
			#Check if the read name is in the error dictionary.
			#If the read name is in the error dict means previous conflict with the
			#read, so add the new read there and do not process.
			if name in p_error:
				p_error[name].append((seq,qual))
				private_counter[3] += 1 #error
			elif name in s_error:
				#print(name.decode("UTF-8"),": error",)
				#OBS! We do not copy s_error therefore error MUST BE APPEND to s_error in the merge
				p_error[name]=[(seq,qual)]
				private_counter[3] += 1 #error
			'''
			#Check if SEQ and QUAL have different lengths
			elif (len(seq) != len(qual)):
				#As in this case it does not mean a discrepance between sequences but a problem in
				#The source file, if the read is present in the shared_data, it will stay and this
				#Read will simply be ignored and not counted at all.
				#print("Error: SEQ and QUAL have different lenghts.")
			else:
			'''
			#Only process if SEQ and QUAL are well-formed
			if (len(seq) == len(qual)):
				#Verify the name starts with @ as a FASTQ read should
				#We need to extract it as a range instead of just the first character to
				#get a byte string and not an int representing the chararcter.
				#if name[0:1]!=b"@":
				#Alternatively we can convert the character into the ASCII int instead.
				assert name[0]==ord("@"), "Something is wrong in " + name.decode("UTF-8")
				'''
				if name[0]!=ord("@"):
					print(name.decode("utf-8"))
					print("Something is wrong")
					exit()
				'''
				#Try to get the read. If not it returns None
				s_read = s_data.get(name, None)
				p_read = p_data.get(name, None)

				#Check if the read name is not in the data dictionary.
				if p_read == None:
					if s_read == None:
						#print(name.decode("UTF-8"),": adding new read",)
						p_data[name]=(seq,qual)
						#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
						private_counter[0] += 1 #added
						#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
					else:
						#print("DEDUPE shared")
						deduped = dedupe(seq, qual, s_read)
						if deduped:
							#Only store the read if it has been extended. In this case, the shared item will be erased.
							#If identical, leave the shared copy and add nothing to p_data.
							if(deduped[0] == 2):
								p_data[name]=deduped[1]
								try:
									s_data.pop(name)
								except Exception:
									print(name, "Has already been erased by another process. Continue")
									pass
							#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
							private_counter[deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
							#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
						else:
							p_error[name]=[(seq,qual), s_read]
							#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
							private_counter[3] += 1 #error
							#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
							#pass
				else:
					#print("DEDUPE private")
					deduped = dedupe(seq, qual, p_read)
					if deduped:
						#Only rewrite the read if it has been extended. If identical, do nothing.
						if(deduped[0] == 2):
							p_data[name]=deduped[1]
						#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
						private_counter[deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
						#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
					else:
						#TODO this should be added to the existing list
						p_error[name]=[(seq,qual), p_read]
						#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
						private_counter[3] += 1 #error
						#print("Line ", inspect.currentframe().f_lineno, " - ", shared_arr[:])
						
			#print("leaving if (len(seq) == len(qual))")
		#When while reaches EOF
		#else:
			#print(os.getpid(),"EOF reached with while")
			#pass
	#print("leaving function")
	#Save in the shared dictionary, for which we have to lock it.
	#Ideally no duplicate items must be found and the saving process should be smooth.
	#In the case there were ducplicates in the source file which ended up in different processes
	#We will dedupe them and add a single entry into s_data
	with lock_manager: #Keep other processes from saving to s_data
		locked.value = True #Make other processes wait before continuing reading s_data
		#TODO: The saving code goes in here
		#1. Create intersection between shared and private
		duplicate_names = s_data.keys() & p_data.keys()
		#2. Dedupe whatever is in the intersection and update in the private
		for name in duplicate_names:
			deduped = dedupe(p_data[name][0], p_data[name][0], s_data[name])
			if deduped:
				#Only rewrite the read if it has been extended. If identical, do nothing.
				if(deduped[0] == 2):
					p_data[name]=deduped[1]
				private_counter[deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
			else:
				#TODO this should be added to the existing list
				p_error[name]=[(seq,qual), p_read]
				private_counter[3] += 1 #error
			private_counter[0] -= 1 #Substract the read from added as it has been added to either identical, extended or error
		#3. Update the shared with the private data and update the counter
		s_data.update(p_data)
		shared_arr[shared_arr_index] += private_counter[0]
		shared_arr[shared_arr_index + 1] += private_counter[1]
		shared_arr[shared_arr_index + 2] += private_counter[2]
		shared_arr[shared_arr_index + 3] += private_counter[3]
		p_data.clear()
		locked.value = False
	return True

def init_process(shared_arr_, locked_):
	global shared_arr, locked
	shared_arr = shared_arr_ # must be inherited, not passed as an argument
	locked = locked_

def gzip_nochunks_byte_u3():
	global arr_index
	global shared_arr
	
	#print(shared_arr[:])
	#print(arr_index)
	start_time = time.time()
	for filename in filenames:
		print("Loading",filename," in RAM using ",cpus," processes")
		with io.BytesIO() as gzfile:
			with multiprocessing.Pool(cpus) as pool:
				multiple_results=[]
				for i in range(cpus+1):
					#print("sending job: ", i)
					multiple_results.append(pool.apply_async(getChunk,args=(filename,i,cpus)))
					#pool.apply_async(getChunk,args=(filename,my_shared_list,i,cpus))
				#print("Closing pool")
				pool.close()
				while True:
					pFinished=0
					for res in multiple_results:
						pFinished += res.ready()
					print("Loaded: ", round(100*pFinished/len(multiple_results)),"%", end="\r")
					if pFinished == len(multiple_results):
						print("Loaded")
						break
				#print("Joining pool")
				pool.join()
			print("Merging file")
			for result in multiple_results:
				chunk=result.get()
				gzfile.seek(chunk[0])
				gzfile.write(chunk[1])
			c_size = gzfile.seek(0,2)
			gzfile.seek(0)
			with gzip.GzipFile(mode='rb', fileobj=gzfile) as file:
				print("Estimating compression ratio", end="\r")
				file.seek(1000000)
				c_ratio=gzfile.tell()/file.tell()
				file.seek(0)
				print("Estimating compression ratio:",round(c_ratio,2))
				print("Approx. decompressed file: ", round(c_size/c_ratio))
				print("Decompressing and processing")
				i=1
				arr_index=0
				buffersize=int((c_size/c_ratio)/(cpus))
				v_block=None
				print("Chunk size: ",buffersize)
				print (time.strftime("%c"))
				with multiprocessing.Pool(processes=cpus, initializer=init_process, initargs=(shared_arr,locked,)) as pool:
					multiple_results=[]
					
					while True:
						chunk = file.read(buffersize)
						if chunk == None:
							print(i, " - EOF - flushing")
							print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
							print("Error?")
							break
						elif len(chunk) < buffersize:
							if v_block == None:
								print(i," - File read in one chunk")
								v_block = chunk
								print(round(100*i/(cpus+1)), "% - block ready - ",file.tell()," sending job: ", i, end="\r")
								#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								#processReads(v_block, data_manager, error_manager, lock_manager)
								multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, arr_index)))
								arr_index += 4
								break
							else:
								#print(i," - Last chunk")
								v_block += chunk
								print(round(100*i/(cpus+1)), "% - block ready - ",file.tell()," sending job: ", i, end="\r")
								#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								
								#processReads(v_block, data_manager, error_manager, lock_manager)
								#multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, counters[-1])))
								multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, arr_index)))
								arr_index += 4
								break
						elif file.tell() == buffersize:
							#print(i," - Very beginning of the file")
							v_block = chunk
						else:
							x0=0
							x=chunk.find(b"\n+", x0)
							#print(x)
							while True:
								if x==-1:
									#no hit means in the middle of a sequence
									#print("still inside a block")
									v_block += chunk
									break
								else:
									#print("\\n+: ",x)
									if chunk.rfind(b"\n@",x0,x)==-1:
										x0=x+1
										x=chunk.find(b"\n+",x0)
									else:
										#print("first \\n@ before \\n+:  ",chunk.rfind(b"\n@",x0,x))
										newblock = chunk.rfind(b"\n@",x0,x)+1
										v_block += chunk[:newblock]
										print(round(100*i/(cpus+1)), "% - block ready - ",file.tell()," sending job: ", i, end="\r")
										#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
										
										#processReads(v_block, data_manager, error_manager, counters[-1])
										multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, arr_index)))
										arr_index += 4
										v_block = chunk[newblock:]
										i+=1
										break
							#print("inner wall exit")
					print(100,"%",end="\r\n")
					print("Closing pool")
					pool.close()
					print("Waiting for results to be ready")
					
					while True:
						time.sleep(1)
						pFinished=0
						for res in multiple_results:
							pFinished += res.ready()
						reads_total = sum(shared_arr[:])
						if(cpus > 1):
							reads_added = sum(operator.itemgetter(*range(0,len(shared_arr),4))(shared_arr))
							reads_identical = sum(operator.itemgetter(*range(1,len(shared_arr),4))(shared_arr))
							reads_extended = sum(operator.itemgetter(*range(2,len(shared_arr),4))(shared_arr))
							reads_error = sum(operator.itemgetter(*range(3,len(shared_arr),4))(shared_arr))
						else:
							reads_added = shared_arr[0]
							reads_identical = shared_arr[1]
							reads_extended = shared_arr[2]
							reads_error = shared_arr[3]
						
						print(pFinished,"/",len(multiple_results),"Tot: ", reads_total,"Added: ", reads_added," Iden: ", reads_identical," Ext: ", reads_extended," Err: ", reads_error, "Speed: ", round(reads_total/(time.time()-start_time)), " reads/s              ", end="\r")
						if pFinished == len(multiple_results):
							print(pFinished,"/",len(multiple_results),"Tot: ", reads_total,"Added: ", reads_added," Iden: ", reads_identical," Ext: ", reads_extended," Err: ", reads_error, "                                                                             ")
							print("DONE")
							print(shared_arr[:])
							break
					print("Joining pool")
					pool.join()
					print("outer wall exit")
		print("LOADED")
	#Write the output to disk
	'''
	with multiprocessing.Pool(3) as pool:
		pool.map(save4, dict(data),100)
		print("Closing pool")
		pool.close()
		print("Joining pool")
		pool.join()
	'''

	print (time.strftime("%c"))
	print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
	print("Saving to disk")
	#print(timeit.Timer(save1).timeit(number=1))
	#print(timeit.Timer(save2).timeit(number=1))
	#print(timeit.Timer(save3).timeit(number=1))
	#data={}
	gc.collect()
'''
def save1():
	with gzip.open('all_1.fastq.gz', 'wb') as gzfile:
		for key, value in data.items():
			gzfile.write(key+b"\n"+value[0]+b"\n+\n"+value[1]+b"\n")
		gzfile.close()

def save2():
	print("Building file")
	with io.BytesIO() as file:
		for key, value in data.items():
			file.write(key+b"\n"+value[0]+b"\n+\n"+value[1]+b"\n")
		print("Writting file")
		with gzip.open('all_2.fastq.gz', 'wb') as gzfile:
			gzfile.write(file.getvalue())
			gzfile.close()
		file.close()
	print("Done")

def save3():
	gzfileStream = io.BytesIO()
	with gzip.GzipFile(mode='wb', fileobj=gzfileStream) as gzfile:
		for key, value in data.items():
			gzfile.write(key+b"\n"+value[0]+b"\n+\n"+value[1]+b"\n")
		gzfile.close()
		with open('all_3.fastq.gz', 'wb') as fh:
			fh.write(gzfileStream.getvalue())
			fh.close()

def save4(iter):
	print("hjgj: ")
'''
ncbu3=timeit.Timer(gzip_nochunks_byte_u3).timeit(number=1)

print("ncbu3=",ncbu3)

print (time.strftime("%c"))

exit()
