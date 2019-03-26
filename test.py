import shutil
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
import csv
import re

def  humanbytes(size, unit = True):
	#2**10 = 1024
	power = 2**10
	n = 0
	Dic_powerN = {0 : 'bytes', 1: 'KB', 2: 'MB', 3: 'GB', 4: 'TB'}
	while size > power:
		size /=  power
		n += 1
	if unit == True:
		return str(round(size,1)) + " " + Dic_powerN[n]
	else:
		return str(round(size,1))

def get_size(obj, seen=None):
	"""Recursively finds size of objects"""
	size = sys.getsizeof(obj)
	if seen is None:
		seen = set()
	obj_id = id(obj)
	if obj_id in seen:
		return 0
	# Important mark as seen *before* entering recursion to gracefully handle
	# self-referential objects
	seen.add(obj_id)
	if isinstance(obj, dict):
		size += sum([get_size(v, seen) for v in obj.values()])
		size += sum([get_size(k, seen) for k in obj.keys()])
	elif hasattr(obj, '__dict__'):
		size += get_size(obj.__dict__, seen)
	elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
		size += sum([get_size(i, seen) for i in obj])
	return size


def sprint(*args, end="\r\n"):
	global messages
	
	finString = ""
	for arg in args:
		finString += str(arg)
	
	messages.put((end, finString))
	
	return True

def count(cpus, messages):
	#TODO: transform this into a singleton
	global readsCounter
	global counterActive
	last_end = "\r\n"
	columns, rows = os.get_terminal_size(0)
	
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
				
				if messages.empty() == False:
					if last_end == "\r":
						print("\033[1A\033[K", end = "")
					else:
						print("\033[K", end = "")
				
					while messages.empty() == False:
						message = messages.get()
						print("\033[K", end = "")
						print(message[1], end = message[0], flush = True)
						last_end = message[0]
					if last_end == "\r":
						print("")
				print("Tot: ", reads_total,"Added: ", reads_added," Iden: ", reads_identical," Ext: ", reads_extended," Err: ", reads_error, "Speed: ", round((reads_total - start_reads)/(time.time()-start_time)), " reads/s", end="\r")
		else:			
			while messages.empty() == False:
				message = messages.get()
				print("\033[K", end = "")
				print(message[1], end = message[0], flush = True)
				last_end = message[0]
	return

def gzCompress(slice_i, slice_f):
	global s_data
	global s_data_keys
	with io.BytesIO() as gzfileStream:
		with gzip.GzipFile(mode='wb', fileobj=gzfileStream) as file:
			for index in range(slice_i, slice_f):
				file.write(s_data_keys[index]+b"\n"+s_data[s_data_keys[index]][0]+b"\n+\n"+s_data[s_data_keys[index]][1]+b"\n")
			file.close()
		sprint("Waiting to write into the file", end="\r")
		with outputLock:
			sprint("Writing into the file", end="\r")
			with open(outputFile, 'ab') as gzfile:
				gzfile.write(gzfileStream.getvalue())
				gzfile.close()
		gzfileStream.close()
	sprint("Releasing memory: gzCompress", end="\r")
	gc.collect()
	sprint("Done: gzCompress", end="\r")
	return 

def init_gzCompress(outputFile_, outputLock_):
	global outputFile, outputLock
	outputLock = outputLock_
	outputFile = outputFile_


def best_overlap(read1,qual1,read2,qual2, minoverlap=25):
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
	for r1_offset in range(len(read1)-minoverlap):
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



def dedupe(seq, qual, stored_read, minoverlap=25):
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
		consensus = best_overlap(seq,qual, stored_read[0], stored_read[1], minoverlap)
		
		#Check if overlap was not found
		if consensus == -1:
			#Try to find overlap:
			#read2
			#----------------
			#          ---------------
			#                    read1
			consensus = best_overlap(stored_read[0], stored_read[1], seq, qual, minoverlap)
		
		#Check if overlap was not found
		if consensus == -1:
			return None
			
		#If overlap was found
		else:
			#Add consensus (contains SEQ and QUAL) to data dict
			#print("consensus found",)
			return (2, consensus)


#def processReads(zByteString):
def processReads(chunk, minoverlap=25):
	'''
	It adds the reads in byteString to the shared dictionary data
	if they are correct/extended/overlapped or to the error
	dictionary if they were not properly processed.
	Returns true if it finished succesfully or exits otherwise
	'''
	global s_data
	global s_error
	procID = None
	
	file = None
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
	#sprint("v_block received")
	#Counter contains an array("added", "identical", "extended", "error")
	if procID == None:
		print("Counter not initialised. Exiting")
		exit()
	
	#print("Allocated procID: ", procID)
	#print(readsCounter[:])
	#time.sleep(1)
	
	#byteString = zlib.decompress(zByteString)
	#del zByteString
	#sprint("v_block decompressed")
	
	#Load the byteString into a byteStream to loop over it
	with io.BytesIO(chunk) as file:
		#Erase the chunk to release mempory
		chunk_length = len(chunk)
		del chunk
		sprint("Releasing memory: ProcessReads", end="\r")
		gc.collect()
		sprint("Done: ProcessReads", end="\r")
		#Make sure we are at the beginning of the stream.
		file.seek(0)
		
		#Execute until we reach the EOF
		while file.tell() < chunk_length: 
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
					deduped = dedupe(seq, qual, p_read, minoverlap)
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
						#sprint("Added error in dedupe private")
						#print("\n>", name.decode("UTF-8"), "\n", read[0].decode("UTF-8"), file=sys.stderr)
						private_counter[3] += 2 #error
						private_counter[0] -= 1 #error
						#print("Line ", inspect.currentframe().f_lineno, " - ", readsCounter[:])
				else:
				#If not in p_data, check if it's in s_data.
					if s_read:
						#print("DEDUPE shared")
						deduped = dedupe(seq, qual, s_read, minoverlap)
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
							#sprint("Added error in dedupe shared")
							#print("\n>", name.decode("UTF-8"), "\n", s_read[0].decode("UTF-8"), "\n>", name.decode("UTF-8"), "\n", seq.decode("UTF-8"), file=sys.stderr)
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
		#gc.collect(0)
	#print("leaving function")
	#Save in the shared dictionary, for which we have to lock it.
	#Ideally no duplicate items must be found and the saving process should be smooth.
	#In the case there were ducplicates in the source file which ended up in different processes
	#We will dedupe them and add a single entry into s_data
	#print("\n", os.getpid(), " - All reads processed.")
	#print(readsCounter[:])
	readsCounter[procID + 1] += private_counter[0]
	readsCounter[procID + 2] += private_counter[1]
	readsCounter[procID + 3] += private_counter[2]
	readsCounter[procID + 4] += private_counter[3]
	
	#print(readsCounter[:])
	#results = zlib.compress(pickle.dumps((p_data, p_error), protocol=4))
		
	updaterQueue.put_nowait(zlib.compress(pickle.dumps((p_data, p_error), protocol=4),1))
	#updaterQueue.put(pickle.dumps((p_data, p_error), protocol=4))
	
	del p_data
	del p_error
	
	with readsCounterLock:
		readsCounter[procID] = 0
	
	sprint("Releasing memory: Processreads", end="\r")
	gc.collect()
	sprint("Done: Processreads", end="\r")
	#print(procID, " Finished")
	
	return True

def init_processReads(readsCounter_, readsCounterLock_, updaterQueue_):
	#print("init_processReads start")
	global readsCounter, readsCounterLock, FASTQ_chunks, updaterQueue
	readsCounter = readsCounter_ # must be inherited, not passed as an argument
	readsCounterLock = readsCounterLock_
	updaterQueue = updaterQueue_
	#print("init_processReads stop") 


class updateShared:
	#print("updateshared class")
	#TODO: transform this into a singleton
	_updater = None
	_running = False
	
	u_data = dict()
	u_error = dict()
	
	def start(self):
		#print("updateshared start")
		if self._updater:
			if(self._updater.isAlive() == True):
				print(threading.currentThread().getName(), " is already running")
				return
		self.u_data = dict()
		self.u_error = dict()
		self._updater = threading.Thread(target=self._updateShared, name="updateShared")
		#self._updater.setDaemon(True)
		self._updater.start()
		#print("updateshared start stop")
		
	def stop(self):
		#print("updateshared stop start")
		sprint("Waiting for the queues to be empty")
		while updaterQueue.empty() == False:
			sprint("Empty: ", updaterQueue.empty(), " Size: ", end="\r")
			time.sleep(0.1)
		#print("Queues empty")
		sprint("Trying to stop updater", end="\r")
		self._running = False
		#self._updater.join()
		while self._updater.isAlive():
			sprint("Trying to stop updater", end="\r")
			time.sleep(0.1)
		sprint("Releasing memory: Updater", end="\r")
		gc.collect()
		sprint("Done: Updater", end="\r")
		sprint("Updater stopped")

	def getData(self):
		return self.u_data.items()
		
	def getErrors(self):
		return self.u_error.items()
	
	def length(self):
		return (len(self.u_data), len(self.u_error))
	
	def size(self):
		return (sys.getsizeof(self.u_data), sys.getsizeof(self.u_error))
	
	def _updateShared(self):
		#print("_updateshared start")
		#TODO: Pass the values "properly"
		self._running = True
		
		global updaterQueue
		
		#print(threading.currentThread().getName(), " Launched")
		
		while self._running == True:
			sprint("U_files left: ", updaterQueue.qsize(), end="\r")
			while updaterQueue.empty() == False:
				sprint("U_files left: ", updaterQueue.qsize(), end="\r")
				#name, read = updaterQueue.get()
				p_data, p_error = pickle.loads(zlib.decompress(updaterQueue.get()))
				#p_data, p_error = pickle.loads(updaterQueue.get())
				
				for e_name, e_reads in p_error.items():
					#0 Add e_name to self.u_error
					self.u_error.setdefault(e_name,list()).extend(e_reads)
					#readsCounter[(cpus*5) + 3] += 1 #error
					
					#1 If a e_name is in self.u_data, transfer it to self.u_error.
					if e_name in self.u_data:
						self.u_error[e_name].append(self.u_data[e_name])
						self.u_data.pop(e_name)
						readsCounter[(cpus*5) + 3] += 1 #error
						readsCounter[(cpus*5)] -= 1 #Substract the read from added
						#sprint("Added error in e_name self.u_data  \n>")
						#for read in self.u_error[e_name]:
						#	sprint(">", e_name, "\n", read[0])
				
				for name, read in p_data.items():
					#2. If name is in self.u_error, transfer read to self.u_error
					if name in self.u_error:
						self.u_error[name].append(read)
						readsCounter[(cpus*5) + 3] += 1 #error
						readsCounter[(cpus*5)] -= 1 #Substract the read from added
						#sprint("Added error in self.u_error name   \n>")
						#for read in self.u_error[name]:
						#	sprint(">", name, "\n", read[0])
						
					#At this stage, names in self.u_data are guaranteed not to be in  self.u_error
					#3. Dedupe whatever is in self.u_data and update the read value
					elif name in self.u_data:
						#sprint("Found duplicate: ", name)
						deduped = dedupe(read[0], read[1], self.u_data[name], minoverlap)
						if deduped:
							#Only rewrite the read if it has been extended.
							if(deduped[0] == 2):
								p_data[name]=deduped[1]
							readsCounter[(cpus*5) + deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
							#private_counter[deduped[0]] += 1 #dedupe() will tell if it is 1 (identical) or 2 (extended)
						else:
							self.u_error.setdefault(name,list()).extend([self.u_data[name], read])
							self.u_data.pop(name)
							readsCounter[(cpus*5) + 3] += 2 #error
							#sprint("Added error in self.u_data:")
							#for read in self.u_error[name]:
							#	print("\n>", name.decode("UTF-8"), "\n", read[0].decode("UTF-8"), file=sys.stderr)
						readsCounter[(cpus*5)] -= 1 #or is 2? #Substract the read from added as it has been added to either identical, extended or error
					else:
						#print(name, read)
						self.u_data[name] = read
				#print("empty")
				sprint("Releasing memory: Updater", end="\r")
				gc.collect()
				sprint("Done: Updater", end="\r")
				sprint("U_files left: ", updaterQueue.qsize(), end="\r")
			sprint("U_files left: ", updaterQueue.qsize(), end="\r")
			time.sleep(1)
		sprint(threading.currentThread().getName(), " Leaving")
		return


def decompressChunks(filepath, buffersize):
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
	#We use ceil to round up the size of the buffer and avoid the creation of cpus+1 chunks
	#As we are estimating the compression ratio with the first 10MB of uncompressed data it
	#can always happen that we underestimate the size and a cpus+1 chunks are generated.
	#If the chunk would be more than 1GB, limit it to 1GB to avoid pickling issues
	#file = io.BytesIO(rawArray)
	with open(filepath, 'rb') as gzfile:
		sprint("Estimating compression ratio")
		#Seek the first 10MB (10485760 bytes) of uncompressed data
		gz_size = gzfile.seek(0, 2)
		gzfile.seek(0)
		
		sprint("Compressed file size: ", humanbytes(gz_size))

		with gzip.open(gzfile, 'rb') as file:
		
			file.seek(10*1024*1024)
			
			gz_ratio = gzfile.tell()/file.tell()
			
			f_size = gz_size/gz_ratio
			
			file.seek(0)
			gzfile.seek(0)
			
			sprint("Estimated compression ratio:", round(gz_ratio, 2))
			sprint("Approx. decompressed file size: ", humanbytes(f_size))
			
			if (f_size/cpus) < buffersize:
				buffersize = math.ceil(f_size/cpus)
				
			sprint("Using chunk size of: ", humanbytes(buffersize))
			
			sprint("Decompressing and generating chunks")
			
			offset_0 = 0
			offset = buffersize
			peek_size = 10*1024
			chunks = list()
			peek_buffer = b''
			with io.BytesIO() as content:
				while True:
					sprint("Start while 1", end="\r")
					if offset_0 == 0 and offset >= f_size:
						content.write(peek_buffer)
						content.write(file.read())
						chunks.append(content.getvalue())
						content.close()
						file.close()
						gzfile.close()
						sprint("Releasing memory: DecompressChunks", end="\r")
						gc.collect()
						sprint("Done: DecompressChunks", end="\r")
						sprint("File read in one chunk")
						return chunks
					else:
						sprint(" + Erasing content", end="\r")
						content.seek(0)
						content.truncate()
						sprint(" + Emptying peek in content", end="\r")
						sprint(" + old peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
						sprint(" + old content: ", content.getvalue()[:11], " ... ", content.getvalue()[-10:], end="\r")
						content.write(peek_buffer[:buffersize])
						sprint(" + new content: ", content.getvalue()[:11], " ... ", content.getvalue()[-10:], end="\r")
						peek_buffer = peek_buffer[buffersize:]
						sprint(" + new peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
						sprint(" + offset_0: ", offset_0, " offset: ", offset, end="\r")
						sprint(" + old tell: ", file.tell(), "/", f_size, end="\r")
						#sprint(" + Reading ", (offset - offset_0) - content.tell(), end="\r")
						sprint(" + Reading ", humanbytes(offset_0), " to ", humanbytes(offset), " / ", humanbytes(f_size), end="\r")
						content.write(file.read((offset - offset_0) - content.tell()))
						sprint(" + new tell: ", file.tell(), "/", f_size, end="\r")
						sprint(" + new content2: ", content.getvalue()[:11], " ... ", content.getvalue()[-10:], end="\r")
						peek_buffer += file.read(peek_size - len(peek_buffer))
						sprint(" + new peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
						sprint(" + new tell2: ", file.tell(), "/", f_size, end="\r")
						seed = peek_buffer.find(b'\n+')+1
						sprint(" + new seed: ", seed, end="\r")
						while f_size > file.tell() and len(peek_buffer) > 0:
							sprint(" + Start while2", end="\r")
							if seed == 0:
								sprint(" + - seed 0", end="\r")
								sprint(" + - - old tell: ", file.tell(), "/", f_size, end="\r")
								sprint(" + - - old peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
								peek_buffer += file.read(peek_size)
								sprint(" + - - new peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
								sprint(" + - - new tell: ", file.tell(), "/", f_size, end="\r")
								seed = peek_buffer.find(b'\n+')+1
								sprint(" + - - new seed: ", seed, end="\r")
							else:
								sprint(" + - seed found: ", seed, end="\r")
								newblock = peek_buffer[:seed].rfind(b'\n@')+1
								sprint(" + - - new block: ", newblock, end="\r")
								if newblock == 0:
									sprint(" + - - newblock 0", end="\r")
									sprint(" + - - - old offset: ", offset, end="\r")
									offset += seed
									sprint(" + - - - new offset: ", offset, end="\r")
									sprint(" + - - - old content: ", content.getvalue()[:11], " ... ", content.getvalue()[-10:], end="\r")
									content.write(peek_buffer[:seed])
									sprint(" + - - - new content: ", content.getvalue()[:11], " ... ", content.getvalue()[-10:], end="\r")
									sprint(" + - - - old peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
									peek_buffer = peek_buffer[seed:]
									sprint(" + - - - new peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
									seed = peek_buffer.find(b'\n+')+1
									sprint(" + - - new seed: ", seed, end="\r")
								else:
									sprint(" + - - newblock found: ", newblock, end="\r")
									sprint(" + - - - old offset: ", offset, end="\r")
									offset += newblock
									sprint(" + - - - new offset: ", offset, end="\r")
									sprint(" + - - - old content: ", content.getvalue()[:11], " ... ", content.getvalue()[-10:], end="\r")
									content.write(peek_buffer[:newblock])
									sprint(" + - - - new content: ", content.getvalue()[:11], " ... ", content.getvalue()[-10:], end="\r")
									sprint(" + - - - old peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
									peek_buffer = peek_buffer[newblock:]
									sprint(" + - - - new peek: ", peek_buffer[:11], " ... ", peek_buffer[-10:], end="\r")
									chunks.append(content.getvalue())
									sprint(" + - - - chunks: ", len(chunks), end="\r")
									sprint(" + - - - old offset_0: ", offset_0, end="\r")
									offset_0 = offset
									sprint(" + - - - new offset_0: ", offset_0, end="\r")
									sprint(" + - - - old offset: ", offset, end="\r")
									offset = offset_0 + buffersize
									sprint(" + - - - new offset: ", offset, end="\r")
									break
						else:
							sprint("EOF reached", end="\r")
							content.write(peek_buffer)
							content.write(file.read())
							chunks.append(content.getvalue())
							content.close()
							file.close()
							gzfile.close()
							sprint("Releasing memory: DecompressChunks", end="\r")
							gc.collect()
							sprint("Done: DecompressChunks", end="\r")
							sprint("chunks: ", len(chunks))
							return chunks
	return False
			
if __name__ == '__main__':
	#from multiprocessing.process import current_process
	#current_process()._config["tempdir"] = "/dev/shm"

	py = psutil.Process(os.getpid())

	parser = argparse.ArgumentParser(description='Dedupe.')
	parser.add_argument('--threads','-t', type=int, default=int(multiprocessing.cpu_count()/2), help='threads to use')
	parser.add_argument('--output','-o', type=str, help='output filename')
	parser.add_argument('--minoverlap', type=int, default=int(25), help='minimum overlap (default: 25)')
	parser.add_argument('--buffersize', type=int, default=int(1024*1024*1024), help='maximum buffer size in bytes')
	parser.add_argument('filepaths', nargs='+', help='files to process')

	args = parser.parse_args()
	cpus = args.threads
	minoverlap = args.minoverlap
	buffersize = args.buffersize
	filepaths = args.filepaths
	
	
	messages = multiprocessing.Queue()
	
	readsCounter = multiprocessing.RawArray(ctypes.c_int,[0]*((cpus*5) + 4))
	readsCounterLock = multiprocessing.Lock()
	
	counterActive = multiprocessing.RawValue(ctypes.c_bool, False)
	counter = multiprocessing.Process(target=count, args=(cpus,messages))
	counter.daemon=True
	counter.start()
	
	updaterQueue = multiprocessing.Queue()
	updater = updateShared()
	
	s_data = dict()
	s_data_keys = None
	s_error = dict()
	
	outputLock = multiprocessing.Lock()
	output = pathlib.Path(args.output)
	output_dir = str(output.parent)
	if(os.path.isdir(output_dir) == False):
		os.mkdir(output_dir)
	output_fileext = ''.join(output.suffixes)
	if "gz" in output_fileext:
		sprint("Gzip file detected")
	elif output_fileext == "":
		sprint("no output file extension, exiting")
		exit()
	output_filename = output.name[:(-len(output_fileext))]
	working_dir = output_dir + "/" + output.name[:(-len(''.join(output.suffixes)))]
	if(os.path.isdir(working_dir) == False):
		os.mkdir(working_dir)
	error_dir = working_dir + "/" + "error_reads"
	if(os.path.isdir(error_dir) == False):
		os.mkdir(error_dir)
	logFile = working_dir + "/" + output_filename + ".log"
	
	sprint (time.strftime("%c"))
	
	previous_filenames = set()
	for filepath in filepaths:
		filename = pathlib.Path(filepath).name
		sprint("Processing file: ", filename)
		#Create a pool of cpus+1 processes to accomodate the extra chunk in case it happens
		#We pass the array to store the reads counted with an initialiser function so the different
		#processes get it by inheritance and not as an argument
		with multiprocessing.Pool(processes=cpus, maxtasksperchild=1, initializer=init_processReads, initargs=(readsCounter, readsCounterLock, updaterQueue)) as pool:
			#Initialise the list that will contain the results of all the processes
			multiple_results=[]
			
			sprint("RAM: ", humanbytes(py.memory_info().rss))
			for chunk in decompressChunks(filepath, buffersize):
				#sprint("Sending job ", len(multiple_results)+1)
				#Send a new process for the chunk
				#sprint(chunk)
				#byteString = memoryview(loadedFASTQ.value)[chunk[0]:chunk[1]] #chunk[1]+1 to include the last character which is \n
				#sprint(byteString[0:10].tobytes(), " -------- ", byteString[-10:-1].tobytes())
				multiple_results.append(pool.apply_async(processReads,args=(chunk, minoverlap)))
				if counterActive.value == False:
					counterActive.value = True
					updater.start()
			sprint(len(multiple_results), " jobs succesfully sent")
			#print("Closing pool")
			pool.close()
			sprint("Releasing memory: Main thread", end="\r")
			gc.collect()
			sprint("Done: Main thread", end="\r")
			sprint("Waiting for results to be ready")
			
			pFinished = [0]*len(multiple_results)
			sprint(sum(pFinished), "/", len(pFinished), end="\r")
			
			while True:
				for index, res in enumerate(multiple_results):
					if((res.ready()==True) and (pFinished[index] == 0)):
						pFinished[index] = int(res.ready())
						res.get() #Force errors in processes to be shown
						sprint(sum(pFinished), "/", len(pFinished), end="\r")
				if sum(pFinished) == len(multiple_results):
					sprint(sum(pFinished), "/", len(pFinished), " DONE")
					#sprint("Joining pool")
					pool.join()
					break
				else:
					time.sleep(1)
			
			#sprint(readsCounter[:])
			del multiple_results
		
		sprint("Releasing memory: Main thread", end="\r")
		gc.collect()
		sprint("Done: Main thread", end="\r")
		updater.stop()
		counterActive.value = False
		
		sprint("Updating the shared dictionaries")
		sprint("Length: ", updater.length())
		sprint("Updating valid data")
		s_data.update(updater.getData())
		sprint("Updating errors")
		for Ename, Ereads in updater.getErrors():
			#If Ename in s_data then it cannot be in s_error.
			
			if Ename in s_data:
				#Transfer the read to s_error
				s_error.setdefault(Ename,dict()).setdefault(s_data[Ename],set()).update(previous_filenames)
				s_data.pop(Ename)
				
				#The Eread(s) in updater.getErrors() needs to be added to s_error too.
				for Eread in Ereads:
					s_error[Ename].setdefault(Eread, set()).add(filename)
				
				#error increases in 1 for the transfered read from s_data. The one in updater.getErrors() was already taken into account
				readsCounter[(cpus*5) + 3] += 1
				
				#Substract the read from added as we transfered it to error
				readsCounter[(cpus*5)] -= 1
			
			elif Ename in s_error:
				#If Ename in s_error we need to check if the read already existed
				for Eread in Ereads:
					deduped = None
					for s_Eread in s_error[Ename].keys():
						deduped = dedupe(Eread[0], Eread[1], s_Eread, minoverlap)
						if deduped:
							#If Eread is either identical or contained in s_Eread
							if(deduped[0] == 1):
								s_error[Ename][s_Eread].add(filename)
							#If Eread and s_Eread have been extended.
							elif(deduped[0] == 2):
								#Add the deduped read to the s_error with all the filenames from the s_Eread
								s_error[Ename].setdefault(deduped[1], set()).update(s_error[Ename][s_Eread])
								#Erased s_Eread as it has been superceeded by the deduped one
								s_error[Ename].pop(s_Eread)
								#Add the filename for Eread to s_error
								s_error[Ename][deduped[1]].add(filename)
							else:
								print("Error!!")
								exit()
							break
					if deduped == None:
						#Eread is totally new for Ename, add it
						s_error[Ename].setdefault(Eread, set()).add(filename)
			else:
				#If Ename not in s_data or s_error, error appeared within the file itself (p_data or u_data)
				#Add a totally new entry to s_error
				for Eread in Ereads:
					s_error.setdefault(Ename,dict()).setdefault(Eread, set()).add(filename)
		sprint("Shared dictionaries updated")
		sprint("File processed succesfully")
		previous_filenames.add(filename)
		sprint("Releasing memory: Main thread", end="\r")
		gc.collect()
		sprint("Done: Main thread", end="\r")
		
	#Write the output to disk
	sprint("S_data: ", len(s_data))
	
	sprint("Exporting the keys")
	
	s_data_keys = multiprocessing.RawArray(ctypes.c_char_p, list(s_data.keys()))
	
	sprint("Compressing the results")
	
	with multiprocessing.Pool(cpus, initializer=init_gzCompress, initargs=(str(output), outputLock)) as pool:
		multiple_results = []
		
		for i in range(cpus):
			multiple_results.append(pool.apply_async(gzCompress,args=((i*len(s_data))//cpus, ((i+1)*len(s_data))//cpus)))
		
		sprint("Closing pool")
		pool.close()

		sprint("Waiting for the compressed results to be ready")
					
		pFinished = [0]*len(multiple_results)
		sprint(sum(pFinished), "/", len(pFinished), end="\r")
		
		while True:
			for index, res in enumerate(multiple_results):
				if((res.ready()==True) and (pFinished[index] == 0)):
					res.get()
					pFinished[index] = int(res.ready())
					sprint(sum(pFinished), "/", len(pFinished), end="\r")
			if sum(pFinished) == len(multiple_results):
				pool.join()
				sprint(sum(pFinished), "/", len(pFinished), " DONE")
				break
			else:
				time.sleep(1)
		#sprint(readsCounter[:])
		del multiple_results
	'''
	s_data_keys = None
	s_data.clear()
	del s_data_keys
	del s_data
	sprint("Releasing memory: Main thread", end="\r")
	gc.collect()
	sprint("Done: Main thread", end="\r")
	'''
	sprint("Generating error table")
	
	with open(error_dir + ".csv", 'w', newline='') as csvfile:
		spamwriter = csv.writer(csvfile, delimiter=';',	quotechar='"', quoting=csv.QUOTE_MINIMAL)
		spamwriter.writerow(["Read name", "Read"]+list(previous_filenames))
		for Ename, Ereads in s_error.items():
			for Eread, Efiles in Ereads.items():
				row = [Ename.decode("UTF-8"), Eread[0].decode("UTF-8")]
				for filename in previous_filenames:
					if filename in Efiles:
						row.append("X")
					else:
						row.append(" ")
				spamwriter.writerow(row)
	sprint("Generating error files")
	error_files = dict()
	for filename in previous_filenames:
		error_files[filename] = gzip.open(error_dir + "/" + filename, 'wb')
	sprint(error_files)
	for Ename, Ereads in s_error.items():
			for Eread, Efiles in Ereads.items():
				for Efile in Efiles:
					error_files[Efile].write(Ename+b"\n"+Eread[0]+b"\n+\n"+Eread[1]+b"\n")
	del s_error	
	
	for filename in error_files:
		error_files[filename].close()
	
	sprint (time.strftime("%c"))
	sprint("RAM: ", humanbytes(py.memory_info().rss))
	


