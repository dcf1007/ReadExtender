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
import math
#from multiprocessing.process import current_process
#current_process()._config["tempdir"] = "/dev/shm"

print (time.strftime("%c"))

py = psutil.Process(os.getpid())

parser = argparse.ArgumentParser(description='Dedupe.')
parser.add_argument('--threads','-t', type=int, default=int(multiprocessing.cpu_count()/2), help='threads to use')
parser.add_argument('filenames', nargs='+', help='files to process')

args = parser.parse_args()
cpus = args.threads
filenames = args.filenames


readsCounter = multiprocessing.RawArray(ctypes.c_int,[0]*(cpus+2)*4)

data_queue = multiprocessing.Queue()
error_queue = multiprocessing.Queue()
s_data = dict()
s_error = dict()

#lock_manager = manager.Lock()
#locked = multiprocessing.RawValue(ctypes.c_bool, False)

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
		#return (2,stored_read)
		return (1,) #As it is contained in the stored read, there is no need to do anything. Considered identical
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


def processReads(byteString, procID = None):
	'''
	It adds the reads in byteString to the shared dictionary data
	if they are correct/extended/overlapped or to the error
	dictionary if they were not properly processed.
	Returns true if it finished succesfully or exits otherwise
	'''
	global s_data
	global s_error
	
	p_error={}
	p_data={}
	
	private_counter = [0,0,0,0]

	#Counter contains an array("added", "identical", "extended", "error")
	if procID == None:
		print("Counter not initialised. Exiting")
		exit()
	
	#Load the byteString into a byteStream to loop over it
	with io.BytesIO(byteString) as file:

		#Make sure we are at the beginning of the stream.
		file.seek(0)
		
		#Execute until we reach the EOF
		while file.tell() < len(byteString): 
			#TODO: Explanation of internal counter. Basically to avoid bottleneck by writting too often to a shared object
			if(sum(private_counter) % 1000 == 0):
				readsCounter[procID] += private_counter[0]
				readsCounter[procID + 1] += private_counter[1]
				readsCounter[procID + 2] += private_counter[2]
				readsCounter[procID + 3] += private_counter[3]
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
			
			#Check if the read name is in the error dictionary.
			#If the read name is in the error dict means previous conflict with the
			#read, so add the new read there and do not process.
			if name in p_error:
				p_error[name].append((seq,qual))
				private_counter[3] += 1 #error
			elif name in s_error:
				#print(name.decode("UTF-8"),": error",)
				#OBS! We do not copy s_error therefore error MUST BE APPEND to s_error in the merge
				p_error[name] = s_error[name]
				p_error[name].append((seq,qual))
				
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
				#s_read = s_data.get(name, None) <--- We are replacing this with the dataNames to accelerate the reading
				p_read = p_data.get(name, None)
				s_read = s_data.get(name, None)
				#Check if the read name is in p_data dictionary.
				if p_read:
					#print("DEDUPE private")
					deduped = dedupe(seq, qual, p_read)
					if deduped:
						#Only rewrite the read if it has been extended. If identical, do nothing.
						if(deduped[0] == 2):
							p_data[name]=deduped[1]
						#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
						private_counter[deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
						#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
					else:
						#TODO this should be added to the existing list
						p_error[name]=[(seq,qual), p_read]
						p_data.pop(name)
						#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
						private_counter[3] += 1 #error
						#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
				else:
				#If not in p_data, check if it's in s_data.
					if s_read:
						#print("DEDUPE shared")
						deduped = dedupe(seq, qual, s_read)
						if deduped:
							#Only store the read if it has been extended. In this case, the shared item will be erased.
							#If identical, leave the shared copy and add nothing to p_data.
							if(deduped[0] == 2):
								p_data[name]=deduped[1]
							#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
							private_counter[deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
							#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
						else:
							p_error[name]=[(seq,qual), s_read]
							#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
							private_counter[3] += 1 #error
							#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
							#pass
					else:
					#If not in p_data or s_data, it is a new read.
						#print(name.decode("UTF-8"),": adding new read",)
						p_data[name]=(seq,qual)
						#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
						private_counter[0] += 1 #added
						#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
						
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
	#print("\n", os.getpid(), " - All reads processed. Waiting signal to add them to the shared dictionary")

	readsCounter[procID] += private_counter[0]
	readsCounter[procID + 1] += private_counter[1]
	readsCounter[procID + 2] += private_counter[2]
	readsCounter[procID + 3] += private_counter[3]
		
	return p_data, p_error

def init_processReads(readsCounter_):
	global readsCounter
	readsCounter = readsCounter_ # must be inherited, not passed as an argument

def gzip_nochunks_byte_u3():
	global readsCounter
	global s_data #Global shared dictionary containing the deduped reads
	global s_error #Global shared dictionary containing the error reads

	start_time = time.time()
	for filename in filenames:
		print("Loading ",filename," in RAM using ",cpus," processes")
		with io.BytesIO() as gzfile:
			with multiprocessing.Pool(cpus+1) as pool:
				multiple_results=[]
				for i in range(cpus+1):
					print("sending job: ", i, end="\r")
					multiple_results.append(pool.apply_async(getChunk,args=(filename,i,cpus)))
					#pool.apply_async(getChunk,args=(filename,my_shared_list,i,cpus))
				print("")
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
				#Seek the first 10MB (10485760 bytes) of uncompressed data
				#file.seek(10485760)
				file.seek(1000000)
				
				#Read the position of the pointers in the compressed and uncompressed data to estimate
				#the compression ratio
				c_ratio=gzfile.tell()/file.tell()
				
				print("Estimating compression ratio:",round(c_ratio,2))
				print("Approx. decompressed file: ", round(c_size/c_ratio), " bytes")
				
				#We use ceil to round up the size of the buffer and avoid the creation of cpus+1 chunks
				#As we are estimating the compression ratio with the first 10MB of uncompressed data it
				#can always happen that we underestimate the size and a cpus+1 chunks are generated.
				buffersize = math.ceil(c_size / (c_ratio * cpus))
				#buffersize = int((c_size/c_ratio)/(cpus))
				print("Using chunk size of: ", buffersize)
				
				
				print("Decompressing and processing")
				
				#Seek the start of the file to start the real decompression
				file.seek(0)
				
				#Initialise v_block
				v_block=None
				
				print (time.strftime("%c"))
				
				#Create a pool of cpus+1 processes to accomodate the extra chunk in case it happens
				#We pass the array to store the reads counted with an initialiser function so the different
				#processes get it by inheritance and not as an argument
				with multiprocessing.Pool(processes=(cpus+1), initializer=init_processReads, initargs=(readsCounter,)) as pool:
					
					#Initialise the list that will contain the results of all the processes
					multiple_results=[]
					
					#Initialise porcessID number to keep track of the offset for each process to write in the shared arrays.
					arr_index=0
					
					#This loop allows to scan the whole file using chunks of defined size
					while True:
						#Read a chunk of length buffersize
						chunk = file.read(buffersize)
						
						#Chunk should be at least an empty byte string b''. Something went wrong if it's None
						if chunk == None:
							#print(int((arr_index/4)+1), " - EOF - flushing")
							#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
							print("Error?")
							exit()
							break
						#If chunk length is smaller than buffer size it's either the last chunk or the whole file read at once 
						elif len(chunk) < buffersize:
							#If no previous v_block was loaded, this is the first and last chunk and therefore the whole file
							if v_block == None:
								
								#Define the chunk as v_block
								v_block = chunk
								
								#Send a new process for the v_block
								multiple_results.append(pool.apply_async(processReads,args=(v_block, arr_index)))
								
								#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								print(int((arr_index/4)+1)," - File read in one chunk")
								print(round(100*int((arr_index/4)+1)/(cpus+1)), "% - block ready - ",file.tell()," sending job: ", int((arr_index/4)+1), end="\r")
								
								#Increase the offset of the shared array for the next process
								arr_index += 4
								
								break
								
							#If there was a previous v_block, this is the last chunk of the file
							else:
								
								#Append the chunk to the existing v_block
								v_block += chunk 
								
								#Send a new process for the v_block
								multiple_results.append(pool.apply_async(processReads,args=(v_block, arr_index)))
								
								#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								print(int((arr_index/4)+1)," - Last chunk")
								print(round(100*int((arr_index/4)+1)/(cpus+1)), "% - block ready - ",file.tell()," sending job: ", int((arr_index/4)+1), end="\r")
								
								#Increase the offset of the shared array for the next process
								arr_index += 4
								
								break
								
						#If the pointer position equals buffersize, we are in the first chunk of a longer file
						elif file.tell() == buffersize:
							#TODO: At this point, v_block should be None. We should assert for this to be sure.
							
							#Define the chunk as v_block
							#At this point, v_block is not ready yet as we don't know if it ends in the middle of a sequence
							v_block = chunk
							
							print(int((arr_index/4)+1)," - Very beginning of the file")
							
						#If it's not at the beginning or the end of the file, we need to define the boundaries of the v_block
						#as the chunk may be ending in the middle of a sequence, and that's not good.
						else:
							'''
							This function works around the structure of a FASTQ read to know whether it is inside a read or not.
							A fastq read is structured in blocks of 4 lines as follows:
								@Name_of_the_read
								SEQUENCE
								+
								QUALITY_VALUES
							Bear in mind that @ and + are valid quality values, so that must be considered when writting the
							algorythms
							
							The first step is find a line that starts with "+" from our offset in the chunk (initially 0).
							That line can either be the line containing + or the line with the quality values if the first
							value is +. That will be our initial seed position.
							
							Then we find the first line starting with "@" contained between the offset position and the seed position.
							The line found is guaranteed to be the beginning of a read. If the quality values of this read would
							start with @, our seed would have stopped at the line that contains the +. So the previous reported line
							starting with @ would still be the beggining of the read.
							
							Once the beggining of the read has been located, that will be the splitting offset to end the previous
							v_block and start a new v_block
							
							'''
							#Set the initial offset to 0
							offset=0
							
							#Find the first line that starts with "+" which will be our seed and store its position
							seed = chunk.find(b"\n+", offset)
							
							#print(seed)
							
							#This loop allows to scan the chunk until the boundaries of the read have been determined
							while True:
								#If no seed was found, we are in the middle of a read and therefore we must continue extending v_block
								if seed==-1:
									#print("still inside a block")
									
									#Append the chunk to the existing v_block
									v_block += chunk
									
									break
								
								#If a seed was found, we need to find the boundaries of the read where the seed is
								else:
									#print("Seed (\\n+) found at: ", seed)
									
									#Find the first line that starts with "@" between the offset and the seed which will be the start of
									#a new read and store its position
									newblock = chunk.rfind(b"\n@", offset, seed)
									
									#If there was no line starting with "@" between the offset and the seed
									if newblock == -1:
										#Set the offset one character after the seed position to start the seed search again
										offset = seed + 1
										
										#Search for a new seed position
										seed = chunk.find(b"\n+", offset)
										
									#If there was a line starting with "@" between the offset and the seed
									else:
										#print("First line starting with @ before the seed: ", newblock)
										
										#Exclude the newline character from the newblock position. At the moment, newblock contains the
										#newline character used in the search, which doesn't correspond to the beginning of the new read
										#and needs to be excluded in the new v_block.
										newblock += 1
										
										#Append the chunk until the newblock position to the existing v_block
										v_block += chunk[:newblock]
										
										#Send a new process for the v_block
										multiple_results.append(pool.apply_async(processReads,args=(v_block, arr_index)))
										
										print(round(100*int((arr_index/4)+1)/(cpus+1)), "% - block ready - ",file.tell()," sending job: ", int((arr_index/4)+1), end="\r")
										#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
										
										#Increase the offset of the shared array for the next process
										arr_index += 4
										
										#Define the rest of the chunk as the new v_block
										v_block = chunk[newblock:]
										
										break
							#print("inner wall exit")
					print(100,"%",end="\r\n")
					print("Closing pool")
					pool.close()
					print("Waiting for results to be ready")
					pFinished = [0]*len(multiple_results)
					while True:
						for index, res in enumerate(multiple_results):
							#if((res.ready()==True) and (pFinished[index] == 0)):
							pFinished[index] = int(res.ready())
						if(cpus > 1):
							reads_added = sum(operator.itemgetter(*range(0,len(readsCounter),4))(readsCounter))
							reads_identical = sum(operator.itemgetter(*range(1,len(readsCounter),4))(readsCounter))
							reads_extended = sum(operator.itemgetter(*range(2,len(readsCounter),4))(readsCounter))
							reads_error = sum(operator.itemgetter(*range(3,len(readsCounter),4))(readsCounter))
						else:
							reads_added = readsCounter[0]
							reads_identical = readsCounter[1]
							reads_extended = readsCounter[2]
							reads_error = readsCounter[3]
						reads_total = sum(readsCounter[:])
						print(sum(pFinished),"/",len(multiple_results),"Tot: ", reads_total,"Added: ", reads_added," Iden: ", reads_identical," Ext: ", reads_extended," Err: ", reads_error, "Speed: ", round(reads_total/(time.time()-start_time)), " reads/s              ", end="\r")
						if sum(pFinished) == len(multiple_results):
							print(sum(pFinished),"/",len(multiple_results),"Tot: ", reads_total,"Added: ", reads_added," Iden: ", reads_identical," Ext: ", reads_extended," Err: ", reads_error, "                                                                             ")
							print("DONE")
							print(readsCounter[:])
							break
						time.sleep(1)
					print("Joining pool")
					pool.join()
					u_data={}
					u_error={}
					for index, res in enumerate(multiple_results):
						p_data, p_error = res.get()
						#TODO: The saving code goes in here
						#1. Create intersection between shared and private
						duplicate_names = u_data.keys() & p_data.keys()
						#2. Dedupe whatever is in the intersection and update in the private
						for name in duplicate_names:
							deduped = dedupe(p_data[name][0], p_data[name][0], u_data[name])
							if deduped:
								#Only rewrite the read if it has been extended. If identical, pop from p_data.
								if(deduped[0] == 2):
									p_data[name]=deduped[1]
								else:
									p_data.pop(name)
								readsCounter[((cpus+1)*4) + deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
								#private_counter[deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
							else:
								#TODO this should be added to the existing list
								p_error[name]=[p_data[name], u_data[name]]
								readsCounter[((cpus+1)*4) + 3] += 1 #error
								#private_counter[3] += 1 #error
							readsCounter[((cpus+1)*4)] -= 1 #Substract the read from added as it has been added to either identical, extended or error
							#private_counter[0] -= 1 #Substract the read from added as it has been added to either identical, extended or error
						
						duplicate_Enames = u_data.keys() & p_error.keys()
						for Ename in duplicate_Enames:
							u_data.pop(name)
							
						#3. Update the shared with the private data and update the counter
						u_data.update(p_data)
						u_error.update(p_error)
						p_data.clear()
						p_error.clear()
						print(index, "Ready", end="\r")
					print("outer wall exit")
					s_data.update(u_data)
					s_error.update(u_error)
		print("LOADED")
	#Write the output to disk
	
	with multiprocessing.Pool(cpus) as pool:
		multiple_results = []
		s_data_keys = list(s_data.keys())
		ks = len(s_data_keys)
		tot = 0
		for i in range(cpus):
			slice = s_data_keys[(i*ks)//cpus:((i+1)*ks)//cpus]
			multiple_results.append(pool.apply_async(gzCompress,args=(slice,)))
			tot += len(slice)
		print(tot)
		print("Closing pool")
		pool.close()
		print("Joining pool")
		pool.join()
		with open('test10b.fastq.gz', 'wb') as fh:
			for res in multiple_results:
				fh.write(res.get())
			fh.close()
	
	print (time.strftime("%c"))
	print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
	#print("Saving to disk")
	#print(timeit.Timer(save1).timeit(number=1))
	#print(timeit.Timer(save2).timeit(number=1))
	#print(timeit.Timer(save3).timeit(number=1))
	#data={}
	gc.collect()

def save1():
	with gzip.open('all_1.fastq.gz', 'wb') as gzfile:
		for key, value in s_data.items():
			gzfile.write(key+b"\n"+value[0]+b"\n+\n"+value[1]+b"\n")
		gzfile.close()

def gzCompress(nameList):
	global s_data
	gzfileStream = io.BytesIO()
	with gzip.GzipFile(mode='wb', fileobj=gzfileStream) as gzfile:
		for name in nameList:
			gzfile.write(name+b"\n"+s_data[name][0]+b"\n+\n"+s_data[name][1]+b"\n")
	return gzfileStream.getvalue()

ncbu3=timeit.Timer(gzip_nochunks_byte_u3).timeit(number=1)

print("ncbu3=",ncbu3)

print (time.strftime("%c"))

exit()
