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

print (time.strftime("%c"))
mypid=os.getpid()
py = psutil.Process(mypid)

parser = argparse.ArgumentParser(description='Dedupe.')
parser.add_argument('--threads','-t', type=int, default=int(multiprocessing.cpu_count()/2), help='threads to use')
parser.add_argument('filenames', nargs='+', help='files to process')
args = parser.parse_args()
cpus = args.threads
filenames = args.filenames

manager = multiprocessing.Manager()
data_manager = manager.dict()
error_manager = manager.dict()
counter_manager=manager.dict({"total":0, "added":0, "identical":0, "extended":0, "error":0})
lock_manager = manager.Lock()

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
		
		assert i*chunksize==fh.tell(), "The pointer in the file and the start of the chunk relative to sourceFile differ. They should be equal"
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

def processReads(byteString, data, error, lock, counter = None):
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

			#FASTQ code is organised in groups of 4 NON-EMPTY lines.
			#3rd line is useless and ignored.
			#safe_readline() makes sure we do not load empty lines on the way.
			name = safe_readline(file) 
			read = safe_readline(file)
			safe_readline(file)
			qual = safe_readline(file)
			
			#If for some reason the file has empty lines at the end of the file, safe_readline will return b''
			if (name==b'' or read==b'' or qual==b''):
				assert name==read==qual==b'', "Error: read not loaded correctly\nName: "+name.decode("UTF-8")+"\nRead: "+read.decode("UTF-8")+"\n Qual: "+qual.decode("UTF-8")
				#print(os.getpid(),": EOF reached with newlines at the end of the file")
				break
				'''
				if name==read==qual==b'':
					print(os.getpid(),": EOF reached with newlines at the end of the file")
				else:
					print("Error: read not loaded correctly")
					exit()
				'''
			#Only process if READ and QUAL are well-formed
			if (len(read) == len(qual)): 
				
				#Lock to prevent different processes writting simultaneously to the dictionary
				with lock:
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
					#Check if the read name is in the error dictionary.
					#If the read name is in the error dict means previous conflict with the
					#read, so add the new read there and do not process.
					#elif name in error:
					if name in error:
						#print(name.decode("UTF-8"),": error",)
						error[name].append((read,qual))
						if counter:
							counter["total"] += 1
							counter["error"] += 1
					#Check if the read name is in the data dictionary. If so, process
					elif name in data:
						#Load the READ and QUAL values stored
						s_read = data[name][0]
						s_qual = data[name][1]
						
						#Check if reads are identical. Omit if so.
						if read == s_read:
							#print(name.decode("UTF-8"),": identical",)
							if counter:
								counter["total"] += 1
								counter["identical"] += 1
							pass
						
						#Check if the new read is contained in the stored read. Omit if so.
						elif read in s_read:
							#print(name.decode("UTF-8"),": s_read is longer")
							if counter:
								counter["total"] += 1
								counter["extended"] += 1
							pass
						
						#Check if the stored read is contained in the new read.
						elif s_read in read:
							#print(name.decode("UTF-8"),": read is longer, storing")
							data[name]=(read,qual)
							if counter:
								counter["total"] += 1
								counter["extended"] += 1
						
						#If reads were neither identical nor contained, try to find a perfect overlap
						else:
							#Try to find overlap:
							#read1
							#----------------
							#          ---------------
							#          read2
							consensus = best_overlap(read,qual,s_read,s_qual)
							
							#Check if overlap was found
							if consensus == -1:
								#Try to find overlap:
								#read2
								#----------------
								#          ---------------
								#          read1
								consensus=best_overlap(s_read,s_qual,read,qual)
							
							#Check if overlap was found
							if consensus == -1:
								
								#If overlap was not found, transfer read to error dict and
								#make sure is erased from data dict too
								if name in error:
									error[name].append((read,qual))
								else:
									error[name]=[(s_read,s_qual),(read,qual)]
									data.pop(name)
								#print(name.decode("UTF-8"),": error",)
								if counter:
									counter["total"] += 1
									counter["error"] += 1
							
							#If overlap was found
							else:
								#Add consensus (contains READ and QUAL) to data dict
								data[name]=consensus
								#print(name.decode("UTF-8"),": consensus found",)
								if counter:
									counter["total"] += 1
									counter["extended"] += 1
					
					#The read was not found in any dictionary. Add to data dict
					else:
						#print(name.decode("UTF-8"),": adding new read",)
						data[name]=(read,qual)
						if counter:
							counter["total"] += 1
							counter["added"] += 1
			
			#If READ and QUAL have different lengths
			else:
				#Transfer read to error dict and make sure is erased from data dict too
				if name in error:
					error[name].append((read,qual))
				else:
					error[name]=[(s_read,s_qual),(read,qual)]
					data.pop(name)
				print("Error: SEQ and QUAL have different lenghts.")
				if counter:
					counter["total"] += 1
					counter["error"] += 1
			#print("leaving if (len(read) == len(qual))")
		#print("leaving while")
		
		#When while reaches EOF
		else:
			#print(os.getpid(),"EOF reached with while")
			pass
	
	#print("leaving function")
	return True


def gzip_nochunks_byte_u3():
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
					pReads=counter_manager["total"]
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
				buffersize=int((c_size/c_ratio)/(cpus))
				v_block=None
				print("Chunk size: ",buffersize)
				print (time.strftime("%c"))
				with multiprocessing.Pool(cpus) as pool:
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
								#print(i," - File read in one chunk")
								v_block = chunk
								print(round(100*i/(cpus+1)), "% - block ready - ",file.tell(),end="\r")
								#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								#print("sending job: ", i)
								#processReads(v_block, data_manager, error_manager, lock_manager)
								multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, lock_manager, counter_manager)))
								break
							else:
								#print(i," - Last chunk")
								v_block += chunk
								print(round(100*i/(cpus+1)), "% - block ready - ",file.tell(),end="\r")
								#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								#print("sending job: ", i)
								#processReads(v_block, data_manager, error_manager, lock_manager)
								multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, lock_manager, counter_manager)))
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
										print(round(100*i/(cpus+1)), "% - block ready - ",file.tell(),end="\r")
										#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
										#print("sending job: ", i)
										#processReads(v_block, data_manager, error_manager, lock_manager)
										multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, lock_manager, counter_manager)))
										v_block = chunk[newblock:]
										i+=1
										break
							#print("inner wall exit")
					print(100,"%",end="\r")
					#print("Closing pool")
					pool.close()
					print("Waiting for results to be ready")
					while True:
						pFinished=0
						for res in multiple_results:
							pFinished += res.ready()
						pReads=counter_manager["total"]
						time.sleep(1)
						print(pFinished,"/",len(multiple_results),"Tot: ", counter_manager["total"],"Added: ", counter_manager["added"]," Iden: ", counter_manager["identical"]," Ext: ", counter_manager["extended"]," Err: ", counter_manager["error"], "Speed: ", counter_manager["total"]-pReads, " reads/s", end="\r")
						if pFinished == len(multiple_results):
							break
					#print("Joining pool")
					pool.join()
					#print("outer wall exit")
		print("LOADED")
		print(len(data_manager))
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
	exit()
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
