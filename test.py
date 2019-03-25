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
parser.add_argument('filenames', nargs='+', help='files to process')
args = parser.parse_args()
filenames=args.filenames

manager = multiprocessing.Manager()
data_manager = manager.dict()
error_manager = manager.dict()
lock_manager = manager.Lock()

i_total=manager.Value("i",0)
i_added=manager.Value("i",0)
i_identical=manager.Value("i",0)
i_extended=manager.Value("i",0)
i_lock=manager.Lock()

def best_overlap(read1,qual1,read2,qual2):
	consensus=-1
	for offset in range(len(read1)):
		if read1[offset:] in read2:
			read1_s=read1.partition(read1[offset:])
			read2_s=read2.partition(read1[offset:])
			if len(read1)-len(read1_s[0])==len(read2)-len(read2_s[2]):
				#offset == len(read1[:offset]))#beginning of read1
				#len(read1)-offset == len(read1[offset:]))#consensus from the read1
				#len(read1)-offset == len(read2[:len(read1[offset:])]))#consensus from the read2
				#len(read2)-len(read1)+offset==len(read2[len(read1[offset:]):]))#end of read2
				consensus = (read1+read2[len(read1)-offset:], qual1+qual2[len(read1)-offset:])
				#print(read_con)
				#print(read1.ljust(len(read_con),b"-"))
				#print(read2.rjust(len(read_con),b"-"))
				break
			else:
				consensus=-1
	return consensus

def gzip_nochunks_byte_u():
	print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
	for filename in filenames:
		print(time.strftime("%c"))
		print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
		print("Loading file: ", filename)
		with gzip.open(filename, 'rb') as fh:
			file=fh.read()
			fh.close()
		fileStream=io.BytesIO(file)
		print (time.strftime("%c"))
		print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
		print("File loaded. Processing")
		while True:
			name = fileStream.readline().strip()
			read = fileStream.readline().strip()
			fileStream.readline()
			qual = fileStream.readline().strip()
			if not (name or read or qual):
				print (time.strftime("%c"))
				print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
				print("EOF reached")
				break
			if (len(read) == len(qual)):
				if name in data:
					s_read = data[name][0]
					s_qual = data[name][1]
					if read == s_read:
						#print(name.decode("UTF-8"),": identical",)
						pass
					elif read in s_read:
						#print(name.decode("UTF-8"),": s_read is longer")
						pass
					elif s_read in read:
						#print(name.decode("UTF-8"),": read is longer, storing")
						data[name]=(read,qual)
					else:
						consensus=best_overlap(read,qual,s_read,s_qual)
						if consensus == -1:
							consensus=best_overlap(s_read,s_qual,read,qual)
						if consensus == -1:
							print(name.decode("UTF-8"),": error",)
						else:
						#	print(name.decode("UTF-8"),": consensus found",)
							data[name]=(consensus)
				else:
					data[name]=(read,qual)
		del file
		del fileStream
		gc.collect()
	#Write the output to disk
	print (time.strftime("%c"))
	print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
	print("Saving to disk")
	#print(timeit.Timer(save1).timeit(number=1))
	print(timeit.Timer(save2).timeit(number=1))
	#print(timeit.Timer(save3).timeit(number=1))
	#data={}
	gc.collect()
	#Done

def gzip_nochunks_byte_u2():
	print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
	for filename in filenames:
		print(time.strftime("%c"))
		print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
		print("Loading file: ", filename)
		with open(filename, "rb") as fh:
			gzfile = fh.read()
			fh.close()
			print (time.strftime("%c"))
			print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
			print("File loaded. Processing")
			gzfileStream=io.BytesIO(gzfile)
			file = gzip.GzipFile(mode="rb", fileobj=gzfileStream)
			while True:
				name = file.readline().strip()
				read = file.readline().strip()
				file.readline()
				qual = file.readline().strip()
				if not (name or read or qual):
					print (time.strftime("%c"))
					print("EOF reached")
					break
				if (len(read) == len(qual)):
					if name in error:
						error.name.append((read,qual))
					elif name in data:
						s_read = data[name][0]
						s_qual = data[name][1]
						if read == s_read:
							#print(name.decode("UTF-8"),": identical",)
							pass
						elif read in s_read:
							#print(name.decode("UTF-8"),": s_read is longer")
							pass
						elif s_read in read:
							#print(name.decode("UTF-8"),": read is longer, storing")
							data[name]=(read,qual)
						else:
							consensus=best_overlap(read,qual,s_read,s_qual)
							if consensus == -1:
								consensus=best_overlap(s_read,s_qual,read,qual)
							if consensus == -1:
								print(name.decode("UTF-8"),": error")
								if error[name]:
									error.name.append((read,qual))
								else:
									error[name]=[(s_read,s_qual),(read,qual)]
									data.pop(name)
							else:
							#	print(name.decode("UTF-8"),": consensus found")
								data[name]=(consensus)
					else:
						data[name]=(read,qual)
			del file
			del gzfileStream
			
		gc.collect()
	#Write the output to disk
	print (time.strftime("%c"))
	print("RAM: ",round(py.memory_info().rss/1024/1024/1024,2))
	print("Saving to disk")
	#print(timeit.Timer(save1).timeit(number=1))
	print(timeit.Timer(save2).timeit(number=1))
	#print(timeit.Timer(save3).timeit(number=1))

	#data={}
	gc.collect()
	#Done

def getChunk(sourceFile, i, cpus):
	with open(sourceFile, "rb") as fh:
		size=os.fstat(fh.fileno()).st_size
		chunksize=round(size/cpus)
		fh.seek(i*chunksize)
		print(i*chunksize==fh.tell())
		return (i*chunksize, fh.read(chunksize))

def safe_readline(byteStream):
	line=b''
	#while (line == b'' and byteStream.tell() < byteStream.getbuffer().nbytes):
	#while (line == b'' and byteStream.tell() < len(byteStream.getbuffer())):
	while (line == b'' and byteStream.tell() < len(byteStream.getvalue())):
		line=byteStream.readline().strip()
	#print(line)
	return line


def processReads(byteString, data, error, lock):
	with io.BytesIO(byteString) as file:
		file.seek(0)
		while file.tell() < len(byteString):
			'''
			if file.tell() >= len(byteStream):
				print(os.getpid(),": EOF by size")
			'''
			name = safe_readline(file)
			read = safe_readline(file)
			safe_readline(file)
			qual = safe_readline(file)
			'''
			if not (name or read or qual):
				if name==read==qual==b'':
					print(os.getpid(),": EOF reached")
					print("sequences: ",len(data))
					break
				else:
					print("error")
					exit()
			'''
			if (name==b'' or read==b'' or qual==b''):
				if name==read==qual==b'':
					print(os.getpid(),": EOF reached with newlines at the end of the file")
				else:
					print("Error: read not loaded correctly")
					exit()

			if (len(read) == len(qual)):
				with lock:
					#We need to extract it as a range instead of just the first character to
					#get a byte string and not an int representing the chararcter.
					#if name[0:1]!=b"@":
					#Alternatively we can convert the character into the ASCII int instead.
					if name[0]!=ord("@"):
						print(name.decode("utf-8"))
						print("Something is wrong")
						exit()
					if name in error:
						#print(name.decode("UTF-8"),": error",)
						error[name].append((read,qual))
					elif name in data:
						s_read = data[name][0]
						s_qual = data[name][1]
						if read == s_read:
							#print(name.decode("UTF-8"),": identical",)
							pass
						elif read in s_read:
							#print(name.decode("UTF-8"),": s_read is longer")
							pass
						elif s_read in read:
						#if s_read in read:
							#print(read)
							#print(s_read)
							#print(name.decode("UTF-8"),": read is longer, storing")
							data[name]=(read,qual)
						else:
						#elif (read != s_read) and (read not in s_read):
							consensus=best_overlap(read,qual,s_read,s_qual)
							if consensus == -1:
								consensus=best_overlap(s_read,s_qual,read,qual)
							if consensus == -1:
								#print(name.decode("UTF-8"),": error",)
								if name in error:
									error[name].append((read,qual))
								else:
									error[name]=[(s_read,s_qual),(read,qual)]
									data.pop(name)
							else:
								#print(name.decode("UTF-8"),": consensus found",)
								data[name]=(consensus)
					else:
						#print(name.decode("UTF-8"),": adding new read",)
						data[name]=(read,qual)
			else:
				print("Error: SEQ and QUAL have different lenghts")
				exit()
			#print("leaving if (len(read) == len(qual))")
		#print("leaving while")
		else: #This else belongs to the while
			print(os.getpid(),"EOF reached with while")
	#print("leaving function")
	return True


def gzip_nochunks_byte_u3():
	for filename in filenames:
		cpus=30
		print("Splitting in ",cpus," cpus")
		with io.BytesIO() as gzfile:
			with multiprocessing.Pool(cpus) as pool:
				multiple_results=[]
				for i in range(cpus+1):
					print("sending job: ", i)
					multiple_results.append(pool.apply_async(getChunk,args=(filename,i,cpus)))
					#pool.apply_async(getChunk,args=(filename,my_shared_list,i,cpus))
				print("Closing pool")
				pool.close()
				print("Joining pool")
				pool.join()
			print("Building file")
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
				print("Decompressing")
				
				i=1
				buffersize=int((c_size/c_ratio)/(cpus))
				v_block=None
				print("Chunk size: ",buffersize)
				
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
								print(i," - File read in one chunk")
								v_block = chunk
								print(i, " - block ready - ",file.tell())
								print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								print("sending job: ", i)
								#processReads(v_block, data_manager, error_manager, lock_manager)
								multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, lock_manager)))
								break
							else:
								print(i," - Last chunk")
								v_block += chunk
								print(i, " - block ready - ",file.tell())
								print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
								print("sending job: ", i)
								#processReads(v_block, data_manager, error_manager, lock_manager)
								multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, lock_manager)))
								break
						elif file.tell() == buffersize:
							print(i," - Very beginning of the file")
							v_block = chunk
						else:
							x0=0
							x=chunk.find(b"\n+", x0)
							#print(x)
							while True:
								if x==-1:
									#no hit means in the middle of a sequence
									print("still inside a block")
									v_block += chunk
									break
								else:
									#print("\\n+: ",x)
									if chunk.rfind(b"\n@",x0,x)==-1:
										x0=x+1
										x=chunk.find(b"\n+",x0)
									else:
										print("first \\n@ before \\n+:  ",chunk.rfind(b"\n@",x0,x))
										newblock = chunk.rfind(b"\n@",x0,x)+1
										v_block += chunk[:newblock]
										print(i, " - block ready - ",file.tell())
										print("----\n",v_block[:10],"...",v_block[-10:],"\n----")
										print("sending job: ", i)
										#processReads(v_block, data_manager, error_manager, lock_manager)
										multiple_results.append(pool.apply_async(processReads,args=(v_block, data_manager, error_manager, lock_manager)))
										v_block = chunk[newblock:]
										i+=1
										break
							print("inner wall exit")
					for res in multiple_results:
						res.wait()
					print("outer wall exit")
		print("Closing pool")
		pool.close()
		print("Joining pool")
		pool.join()
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
	print("hjgj: ",len(iter))

ncbu3=timeit.Timer(gzip_nochunks_byte_u3).timeit(number=1)
ncbu=timeit.Timer(gzip_nochunks_byte_u).timeit(number=1)
ncbu2=timeit.Timer(gzip_nochunks_byte_u2).timeit(number=1)
#ch=timeit.Timer(gzip_chunks).timeit(number=1)

#print("ncb=",ncb)
print("ncbu=",ncbu)
print("ncbu2=",ncbu2)
print("ncbu3=",ncbu3)
#print("ch=",ch)

print (time.strftime("%c"))
#input()
