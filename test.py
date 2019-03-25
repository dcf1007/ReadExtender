import zlib
import gzip 
import time
import timeit
import os
import mmap
import io



print (time.strftime("%c"))

filename = "R2_extended_k-75_p-0.fastq.gz"

def gzip_chunks():
	with gzip.open(filename, 'r') as fh:
		v_block=None
		buffersize=int(27162832335/64)
		i=1
		while True:
			chunk=fh.read(buffersize)
			if not chunk:
				print(i, " - EOF - flushing")
				#print("----\n",v_block,"\n----")
				break
			elif len(chunk) < buffersize:
				print(i, " - File read in one chunk or last chunk")
				v_block = chunk
				#print("----\n",v_block,"\n----")
				break
			elif fh.tell() == buffersize:
				print(i, " - Very beginning of the file")
				v_block = chunk
			elif chunk[:1] == b"@" and v_block[-1:] == b"\n":
					print(i, " - block ready - ",fh.tell())
					#print("----\n",v_block,"\n----")
					v_block = chunk
			else:
				newblock = chunk.find(b"\n@")
				#print("New block: ", newblock)
				if newblock == -1:
					v_block += chunk
					#print("still inside a block")
				else:
					v_block += chunk[:newblock+1]
					print(i, " - block ready - ",fh.tell())
					#print("----\n",v_block,"\n----")
					v_block = chunk[newblock+1:]
			# read line
			'''
			name = fh.readline().strip()
			read = fh.readline().strip()
			fh.readline()
			qual = fh.readline().strip()
			if not (name or read or qual):
				break
			if (len(read) == len(qual)):
				zname=zlib.compress(name)
				# in python 2, print line
				# in python 3
				#print(len(line.strip().encode('utf8')), " vs ", len(zlib.compress(line.strip().encode('utf8'),9)))
				if zname in compressed_data:
					s_read = zlib.decompress(compressed_data[zname][0])
					if read == s_read:
						print("identical")
					elif read in s_read:
						print("read is longer")
					elif s_read in read:
						print("s_read is longer")
					else:
						print("error")
						print(zlib.decompress(zname))
						print(s_read)
						print(read)
				else:
					compressed_data[zname]=(zlib.compress(read),zlib.compress(qual))
					#print(name,read.strip(),qual.strip())
				# check if line is not empty
			'''
		fh.close()

def gzip_chunks2():
	with gzip.open(filename, 'r') as fh:
		file=fh.read()
		fileStream=io.BytesIO(file)
		v_block=None
		buffersize=int(len(file)/64)
		print(len(file))
		i=1
		while True:
			chunk=fileStream.read(buffersize)
			if chunk == None:
				print(i, " - EOF - flushing")
				#print("----\n",v_block,"\n----")
				break
			elif len(chunk) < buffersize:
				print(i," - File read in one chunk or last chunk")
				v_block = chunk
				#print("----\n",v_block,"\n----")
				break
			elif fileStream.tell() == buffersize:
				print(i," - Very beginning of the file")
				v_block = chunk
			elif chunk[:1] == b"@" and v_block[-1:] == b"\n":
					print(i," - block ready - ",fileStream.tell())
					i+=1
					#print("----\n",v_block,"\n----")
					v_block = chunk
			else:
				newblock = chunk.find(b"\n@")
				#print("New block: ", newblock)
				if newblock == -1:
					#print("still inside a block")
					v_block += chunk
				else:
					v_block += chunk[:newblock+1]
					print(i, " - block ready - ",fileStream.tell())
					i+=1
					#print("----\n",v_block,"\n----")
					v_block = chunk[newblock+1:]
			# read line
			'''
			name = fh.readline().strip()
			read = fh.readline().strip()
			fh.readline()
			qual = fh.readline().strip()
			if not (name or read or qual):
				break
			if (len(read) == len(qual)):
				zname=zlib.compress(name)
				# in python 2, print line
				# in python 3
				#print(len(line.strip().encode('utf8')), " vs ", len(zlib.compress(line.strip().encode('utf8'),9)))
				if zname in compressed_data:
					s_read = zlib.decompress(compressed_data[zname][0])
					if read == s_read:
						print("identical")
					elif read in s_read:
						print("read is longer")
					elif s_read in read:
						print("s_read is longer")
					else:
						print("error")
						print(zlib.decompress(zname))
						print(s_read)
						print(read)
				else:
					compressed_data[zname]=(zlib.compress(read),zlib.compress(qual))
					#print(name,read.strip(),qual.strip())
				# check if line is not empty
			'''
		fh.close()


def gzip_nochunks():
	with gzip.open(filename, 'r') as fh:
		file=fh.read()
		fileStream=io.BytesIO(file)
		buffersize=int(len(file)/64)
		print(len(file))

def gzip_nochunks2():
	with open("R2_extended_k-75_p-0.fastq.gz", "r+b") as f:
		mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
		gzfile = gzip.GzipFile(mode="r", fileobj=mm)
		file=gzfile.read()
		fileStream=io.BytesIO(file)
		buffersize=int(len(file)/64)
		print(len(file))

print(timeit.Timer(gzip_nochunks2).timeit(number=1))

nc=timeit.Timer(gzip_nochunks).timeit(number=1)
ch2=timeit.Timer(gzip_chunks2).timeit(number=1)
ch=timeit.Timer(gzip_chunks).timeit(number=1)

print("nc=",nc)
print("ch=",ch)
print("ch2=",ch2)

print (time.strftime("%c"))
#input()
