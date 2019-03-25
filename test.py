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
import pickle
import zlib
import sys
import pathlib
import threading

def sprint(*args, end="\r\n"):
	global messages
	
	columns, rows = os.get_terminal_size(0)
	
	finString = ""
	for arg in args:
		finString += str(arg)
	messages.put(finString.ljust(columns, " ") + end)
	return True

def count(cpus, messages):
	global readsCounter
	global counterActive
	
	while True:
		
		if counterActive.value == True:
			
			start_time = time.time()
			
			if(cpus > 1):
				start_reads = sum(readsCounter[:]) - sum(operator.itemgetter(*range(0,len(readsCounter)-4,5))(readsCounter))
			else:
				start_reads = sum(readsCounter[:]) - readsCounter[0]

			while counterActive.value == True:
				if(cpus > 1):
					reads_added = sum(operator.itemgetter(*range(1,len(readsCounter)-4,5))(readsCounter)) + readsCounter[len(readsCounter)-4]
					reads_identical = sum(operator.itemgetter(*range(2,len(readsCounter)-4,5))(readsCounter)) + readsCounter[len(readsCounter)-3]
					reads_extended = sum(operator.itemgetter(*range(3,len(readsCounter)-4,5))(readsCounter)) + readsCounter[len(readsCounter)-2]
					reads_error = sum(operator.itemgetter(*range(4,len(readsCounter)-4,5))(readsCounter)) + readsCounter[len(readsCounter)-1]
					reads_total = sum(readsCounter[:]) - sum(operator.itemgetter(*range(0,len(readsCounter)-4,5))(readsCounter))
				else:
					reads_added = readsCounter[1]
					reads_identical = readsCounter[2]
					reads_extended = readsCounter[3]
					reads_error = readsCounter[4]
					reads_total = sum(readsCounter[:]) - readsCounter[0]

				while not messages.empty():
					print(messages.get(), end="")
				
				print("Tot: ", reads_total,"Added: ", reads_added," Iden: ", reads_identical," Ext: ", reads_extended," Err: ", reads_error, "Speed: ", round((reads_total - start_reads)/(time.time()-start_time)), " reads/s              ", end="\r")
				
				time.sleep(0.1)
		
		while not messages.empty():
			print(messages.get(), end="")
			
		time.sleep(0.1)
	return

def gzCompress(nameList):
	global s_data
	gzfileStream = io.BytesIO()
	with gzip.GzipFile(mode='wb', fileobj=gzfileStream) as gzfile:
		for name in nameList:
			gzfile.write(name+b"\n"+s_data[name][0]+b"\n+\n"+s_data[name][1]+b"\n")
	return gzfileStream.getvalue()

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
			
			#Calculate r2_offset as the last hit of the query in r2 (the rightmost hit).
			r2_offset = read2.rfind(read1[r1_offset:])
			
			#If there is only one match its value should be 0.
			#Otherwise the query is multiple times in r2.
			if (r2_offset != 0) or (read2.count(read1[r1_offset:]) != 1):
				return consensus
			
			#r1_offset has to be equal to the first hit of the query in r1
			#if there is one match only, otherwise it is repeated.
			if (r1_offset != read1.find(read1[r1_offset:])) or (read1.count(read1[r1_offset:]) != 1):
				return consensus

			#As we have taken the leftmost hit in r1 and the rightmost hit in r2,
			#the sequence after the hit in r1 and before the hit in r2 will only
			#be the same in the case a true consensus is found.
			#OBS! Probably redundant?
			if read1[r1_offset:] == read2[:r2_offset + len(read1)-r1_offset]:
				#print("Only one hit of the query, proceed")
				#OBS!: read1[r1_offset:] == len(read1)-r1_offset
				consensus = (read1+read2[len(read1)-r1_offset:], qual1+qual2[len(read1)-r1_offset:])
				#print(read1.ljust(len(consensus[0]),b"-"))
				#print(read2.rjust(len(consensus[0]),b"-"))
				return consensus
				
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
		
		#Check if overlap was not found
		if consensus == -1:
			#Try to find overlap:
			#read2
			#----------------
			#          ---------------
			#                    read1
			consensus=best_overlap(stored_read[0],stored_read[1],seq,qual)
		
		#Check if overlap was not found
		if consensus == -1:
			return None
			
		#If overlap was found
		else:
			#Add consensus (contains SEQ and QUAL) to data dict
			#print("consensus found",)
			return (2,consensus)


def processReads(byteString):
	'''
	It adds the reads in byteString to the shared dictionary data
	if they are correct/extended/overlapped or to the error
	dictionary if they were not properly processed.
	Returns true if it finished succesfully or exits otherwise
	'''
	global s_data
	global s_error
	procID = None
	
	with readsCounterLock:
		while procID == None:
			for counter in range(0,(len(readsCounter)-4),5):
				if readsCounter[counter] == 0:
					readsCounter[counter] = 1
					procID = counter
					break
	
	p_error={}
	p_data={}
	
	private_counter = [0,0,0,0]

	#Counter contains an array("added", "identical", "extended", "error")
	if procID == None:
		print("Counter not initialised. Exiting")
		exit()
	
	#print("Allocated procID: ", procID)
	#print(readsCounter[:])
	#time.sleep(1)
	
	#Load the byteString into a byteStream to loop over it
	with io.BytesIO(byteString) as file:

		#Make sure we are at the beginning of the stream.
		file.seek(0)
		
		#Execute until we reach the EOF
		while file.tell() < len(byteString): 
			#TODO: Explanation of internal counter. Basically to avoid bottleneck by writting too often to a shared object
			if(sum(private_counter) % 1000 == 0):
				readsCounter[procID + 1] += private_counter[0]
				readsCounter[procID + 2] += private_counter[1]
				readsCounter[procID + 3] += private_counter[2]
				readsCounter[procID + 4] += private_counter[3]
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
				#print("Added error in initial p_error ", name)
			elif name in s_error:
				#We do not append the errors to s_error directly because we identify in which file they were found once the results are returned
				#p_error[name] = s_error[name]
				p_error[name] = [(seq,qual)]
				private_counter[3] += 1 #error
				#print("Added error in initial s_error ", name)
				
			#Check if SEQ and QUAL have different lengths
			elif (len(seq) != len(qual)):
				#As in this case it does not mean a discrepance between sequences but a problem in
				#The source file, if the read is present in the shared_data, it will stay and this
				#Read will simply be ignored and not counted at all.
				print("\n\n\nError: SEQ and QUAL have different lenghts.\n\n\n")
			
			#Only process if SEQ and QUAL are well-formed
			#if len(seq) == len(qual):
			else:
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
						#print("Added error in dedupe private ", name)
						private_counter[3] += 2 #error
						private_counter[0] -= 1 #error
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
							#We don't add the s_read to the p_error as we need to keep track from which file it came from
							#s_read will be transfered after returning the data to main process
							#p_error[name]=[(seq,qual), s_read]
							p_error[name]=[(seq,qual)]
							#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
							private_counter[3] += 1 #error
							#print("Added error in dedupe shared ", name)
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
		del p_read
		del s_read
	#print("leaving function")
	#Save in the shared dictionary, for which we have to lock it.
	#Ideally no duplicate items must be found and the saving process should be smooth.
	#In the case there were ducplicates in the source file which ended up in different processes
	#We will dedupe them and add a single entry into s_data
	#print("\n", os.getpid(), " - All reads processed. Waiting signal to add them to the shared dictionary")
	#print(readsCounter[:])
	readsCounter[procID + 1] += private_counter[0]
	readsCounter[procID + 2] += private_counter[1]
	readsCounter[procID + 3] += private_counter[2]
	readsCounter[procID + 4] += private_counter[3]
	
	#print(readsCounter[:])
	#results = zlib.compress(pickle.dumps((p_data, p_error), protocol=4))
	
	q_data.put(zlib.compress(pickle.dumps((p_data, p_error), protocol=4),2))
	#q_data.put(pickle.dumps((p_data, p_error), protocol=4))
	
	del p_data
	del p_error
	gc.collect()
	
	with readsCounterLock:
		readsCounter[procID] = 0
	
	#print(procID, " Finished")
	
	return True

def init_processReads(readsCounter_, readsCounterLock_, q_data_, q_error_):
	global readsCounter, readsCounterLock
	readsCounter = readsCounter_ # must be inherited, not passed as an argument
	readsCounterLock = readsCounterLock_
	q_data = q_data_
	q_error = q_error_


class updateShared:
	
	_updater = None
	_running = False
	
	u_data = dict()
	u_error = dict()
	
	def start(self):
		if self._updater:
			if(self._updater.isAlive() == True):
				return
		self.u_data = dict()
		self.u_error = dict()
		self._updater = threading.Thread(target=self._updateShared, name="updateShared")
		self._updater.setDaemon(True)
		self._updater.start()
		
	def stop(self):
		#print("Waiting for the queues to be empty")
		while not (q_error.empty() and q_data.empty()):
			#print(int(q_error.qsize() + q_data.qsize()))
			time.sleep(0.1)
		#print("Queues empty")
		#print("trying to stop")
		self._running = False
		self._updater.join()

	def results(self):
		return (self.u_data, self.u_error)
		
	def _updateShared(self):
		#TODO: Pass the values "properly"
		self._running = True
		
		global s_data, s_error
		global q_data, q_error
		
		print(threading.currentThread().getName(), " Launched")
		
		while self._running == True:
					
			while not q_data.empty():
				#name, read = q_data.get()
				p_data, p_error = pickle.loads(zlib.decompress(q_data.get()))
				#p_data, p_error = pickle.loads(q_data.get())
				
				for e_name, e_reads in p_error.items():
					#0 Add e_name to self.u_error
					self.u_error.setdefault(e_name,list()).extend(e_reads)
					readsCounter[(cpus*5) + 3] += 1 #error
					
					#1 If a e_name is in self.u_data, transfer it to self.u_error.
					if e_name in self.u_data:
						self.u_error[e_name].append(self.u_data[e_name])
						self.u_data.pop(e_name)
						readsCounter[(cpus*5) + 3] += 1 #error
						readsCounter[(cpus*5)] -= 1 #Substract the read from added
						#print("Added error in e_name self.u_data ", e_name)
				
				for name, read in p_data.items():
					#2. If name is in self.u_error, transfer read to self.u_error
					if name in self.u_error:
						self.u_error[name].append(read)
						readsCounter[(cpus*5) + 3] += 1 #error
						readsCounter[(cpus*5)] -= 1 #Substract the read from added
						#print("Added error in self.u_error name ", name)
						
					#At this stage, names in self.u_data are guaranteed not to be in  self.u_error
					#3. Dedupe whatever is in self.u_data and update the read value
					elif name in self.u_data:
						#sprint("Found duplicate: ", name)
						deduped = dedupe(read[0], read[1], self.u_data[name])
						if deduped:
							#Only rewrite the read if it has been extended.
							if(deduped[0] == 2):
								p_data[name]=deduped[1]
							readsCounter[(cpus*5) + deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
							#private_counter[deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
						else:
							self.u_error.setdefault(name,list()).append(self.u_data[name], read)
							self.u_data.pop(name)
							readsCounter[(cpus*5) + 3] += 2 #error
							#print("Added error in self.u_data p_data ", name)
						readsCounter[(cpus*5)] -= 1 #or is 2? #Substract the read from added as it has been added to either identical, extended or error
					else:
						#print(name, read)
						self.u_data[name] = read
			
			#print("empty")
			time.sleep(1)
		
		#print(threading.currentThread().getName(), " Leaving")
		return
			
if __name__ == '__main__':

	#from multiprocessing.process import current_process
	#current_process()._config["tempdir"] = "/dev/shm"

	print (time.strftime("%c"))

	py = psutil.Process(os.getpid())

	parser = argparse.ArgumentParser(description='Dedupe.')
	parser.add_argument('--threads','-t', type=int, default=int(multiprocessing.cpu_count()/2), help='threads to use')
	parser.add_argument('--output','-o', type=str, help='output filename')
	parser.add_argument('filepaths', nargs='+', help='files to process')

	args = parser.parse_args()
	cpus = args.threads
	filepaths = args.filepaths

	messages = multiprocessing.Queue()
	
	readsCounter = multiprocessing.RawArray(ctypes.c_int,[0]*((cpus*5) + 4))
	readsCounterLock = multiprocessing.Lock()
	
	counter = multiprocessing.Process(target=count, args=(cpus,messages))
	counterActive = multiprocessing.RawValue(ctypes.c_bool, False)
	counter.start()

	s_data = dict()
	s_data_keys = list()
	s_error = dict()
	
	q_data = multiprocessing.Queue()
	q_error = multiprocessing.Queue()
	
	#lock_manager = manager.Lock()
	#locked = multiprocessing.RawValue(ctypes.c_bool, False)
	
	previous_filenames = []
	for filepath in filepaths:
		filename = pathlib.Path(filepath).name
		sprint("Loading ",filename," in RAM using ",cpus," processes")
		with io.BytesIO() as gzfile:
			with multiprocessing.Pool(cpus+1) as pool:
				multiple_results=[]
				for i in range(cpus+1):
					sprint("sending job: ", i, end="\r")
					multiple_results.append(pool.apply_async(getChunk,args=(filepath,i,cpus)))
					#pool.apply_async(getChunk,args=(filepath,my_shared_list,i,cpus))
				#print("Closing pool")
				pool.close()
				while True:
					pFinished=0
					for res in multiple_results:
						pFinished += res.ready()
					sprint("Loaded: ", round(100*pFinished/len(multiple_results)),"%", end="\r")
					if pFinished == len(multiple_results):
						sprint("Loaded")
						break
				#print("Joining pool")
				pool.join()
				sprint("Merging file")
				for result in multiple_results:
					chunk=result.get()
					gzfile.seek(chunk[0])
					gzfile.write(chunk[1])
				
				del chunk
				del multiple_results
			
			c_size = gzfile.seek(0,2)
			
			gzfile.seek(0)
			
			with gzip.GzipFile(mode='rb', fileobj=gzfile) as file:
				sprint("Estimating compression ratio", end="\r")
				#Seek the first 10MB (10485760 bytes) of uncompressed data
				file.seek(10485760)
				
				#Read the position of the pointers in the compressed and uncompressed data to estimate
				#the compression ratio
				c_ratio=gzfile.tell()/file.tell()
				
				sprint("Estimating compression ratio:",round(c_ratio,2))
				sprint("Approx. decompressed file: ", round(c_size/c_ratio), " bytes")
				
				#We use ceil to round up the size of the buffer and avoid the creation of cpus+1 chunks
				#As we are estimating the compression ratio with the first 10MB of uncompressed data it
				#can always happen that we underestimate the size and a cpus+1 chunks are generated.
				#If the chunk would be more than 1GB, limit it to 1GB to avoid pickling issues
				if (c_size / (c_ratio * cpus)) > 1073741824:
					buffersize = 1073741824
				else:
					buffersize = round(c_size / (c_ratio * cpus))
				sprint("Using chunk size of: ", buffersize)
				
				sprint("Decompressing and processing")
				
				#Seek the start of the file to start the real decompression
				file.seek(0)
				
				#Initialise v_block
				v_block=None
				
				sprint (time.strftime("%c"))
				
				
				#Create a pool of cpus+1 processes to accomodate the extra chunk in case it happens
				#We pass the array to store the reads counted with an initialiser function so the different
				#processes get it by inheritance and not as an argument
				with multiprocessing.Pool(processes=cpus, initializer=init_processReads, initargs=(readsCounter,readsCounterLock, q_data, q_error)) as pool:
					
					#Initialise the list that will contain the results of all the processes
					multiple_results=[]
					
					counterActive.value = True
					
					updater = updateShared()
					updater.start()
					
					#This loop allows to scan the whole file using chunks of defined size
					while True:
						#Read a chunk of length buffersize
						chunk = file.read(buffersize)
						
						#Chunk should be at least an empty byte string b''. Something went wrong if it's None
						if chunk == None:
							#print(int(procID+1), " - EOF - flushing")
							#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
							sprint("Error?")
							exit()
							break
						#If chunk length is smaller than buffer size it's either the last chunk or the whole file read at once 
						elif len(chunk) < buffersize:
							#If no previous v_block was loaded, this is the first and last chunk and therefore the whole file
							if v_block == None:
								
								#Define the chunk as v_block
								v_block = chunk
								
								#Send a new process for the v_block
								multiple_results.append(pool.apply_async(processReads,args=(v_block,)))
								
								#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								sprint( " - File read in one chunk")
								sprint("% - block ready - ", file.tell(), " sending job: ")
								
								break
								
							#If there was a previous v_block, this is the last chunk of the file
							else:
								
								#Append the chunk to the existing v_block
								v_block += chunk 
								
								#Send a new process for the v_block
								multiple_results.append(pool.apply_async(processReads,args=(v_block,)))
								
								#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								
								sprint("% - block ready - ", file.tell(), " sending job: ")
								sprint(" - End of the file")
								
								break
								
						#If the pointer position equals buffersize, we are in the first chunk of a longer file
						elif file.tell() == buffersize:
							#TODO: At this point, v_block should be None. We should assert for this to be sure.
							
							#Define the chunk as v_block
							#At this point, v_block is not ready yet as we don't know if it ends in the middle of a sequence
							v_block = chunk
							
							sprint(" - Very beginning of the file")
							
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
										multiple_results.append(pool.apply_async(processReads,args=(v_block,)))
										
										sprint("% - block ready - ", file.tell(), " sending job: ")
										#print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
										
										#Define the rest of the chunk as the new v_block
										v_block = chunk[newblock:]
										
										break
							#print("inner wall exit")
					sprint("100%")
					del v_block
					del chunk
					file.close()
					gzfile.close()
					
					#print("Closing pool")
					pool.close()
					sprint("Waiting for results to be ready")
					pFinished = [0]*len(multiple_results)
					
					while True:
						for index, res in enumerate(multiple_results):
							if((res.ready()==True) and (pFinished[index] == 0)):
								pFinished[index] = int(res.ready())
								#print(res.get())
								sprint(sum(pFinished), "/", len(multiple_results), end="\r")
						if sum(pFinished) == len(multiple_results):
							sprint("DONE")
							break
						else:
							time.sleep(1)
					sprint("Joining pool")
					pool.join()
					#print("outer wall exit")
					sprint(readsCounter[:])
					updater.stop()
					counterActive.value = False
					sprint("Updating the shared dictionary")
					sprint("Length: ", len(updater.u_data))
					s_data.update(updater.u_data)
					
					for Ename in updater.u_error:
						#1. Check whether the entry exists already in the dictioniary
						#2. If it doesn't exist define a dictionary for it
						if Ename not in s_error:
							s_error[Ename] = dict()
						if Ename in s_data:
							for previous_filename in previous_filenames:
								s_error[Ename].setdefault(previous_filename,list()).append(s_data[Ename])
							s_data.pop(Ename)
							readsCounter[(cpus*5) + 3] += 1 #error
							readsCounter[(cpus*5)] -= 1 #Substract the read from added
							
						s_error[Ename].setdefault(filename, list()).extend(updater.u_error[Ename])
						#u_error.pop(Ename)
						
					#3. Check whether there's an entry for the filename
					#4. If it doesn't exist create the entry and define an empty list on it
					sprint("updated")
					del multiple_results
					gc.collect()
		#sprint("LOADED")
		previous_filenames.append(filename)
	#Write the output to disk
	sprint("S_data: ", len(s_data))
	with multiprocessing.Pool(cpus) as pool:
		multiple_results = []
		s_data_keys = list(s_data.keys())
		ks = len(s_data_keys)
		sprint("S_data_keys: ", len(s_data_keys))
		tot = 0
		for i in range(cpus):
			slice = s_data_keys[(i*ks)//cpus:((i+1)*ks)//cpus]
			multiple_results.append(pool.apply_async(gzCompress,args=(slice,)))
			tot += len(slice)
			del slice
		sprint(tot)
		sprint("Closing pool")
		pool.close()
		sprint("Joining pool")
		pool.join()
		with open(args.output, 'wb') as fh:
			for res in multiple_results:
				#print(gzip.decompress(res.get()))
				fh.write(res.get())
			fh.close()
		s_data_keys = list()
		
		output = pathlib.Path(args.output)
		error_dir = str(output.parent) + "/" + output.name[:(-len(''.join(output.suffixes)))] + "_error"
		if(os.path.isdir(error_dir) == False):
			os.mkdir(error_dir)
		error_files = dict()
		for filename in previous_filenames:
			error_files[filename] = gzip.open(error_dir + "/" + filename, 'wb')
		for name in s_error:
			for filename in s_error[name]:
				for read in s_error[name][filename]:
					error_files[filename].write(name+b"\n"+read[0]+b"\n+\n"+read[1]+b"\n")
		for filename in error_files:
			error_files[filename].close()
	
	sprint (time.strftime("%c"))
	sprint("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
	#print("Saving to disk")
	#print(timeit.Timer(save1).timeit(number=1))
	#print(timeit.Timer(save2).timeit(number=1))
	#print(timeit.Timer(save3).timeit(number=1))
	#data={}
	gc.collect()

	sprint (time.strftime("%c"))
