import zlib
import gzip 
import time
import timeit
import os
import mmap
import io
import gc
import argparse

print (time.strftime("%c"))

parser = argparse.ArgumentParser(description='Dedupe.')
parser.add_argument('filenames', nargs='+', help='files to process')
args = parser.parse_args()
filenames=args.filenames
print(filenames)

def gzip_chunks():
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

compressed_data = {}

def gzip_nochunks_byte():
	with gzip.open(filename, 'rb') as fh:
		file=fh.read()
		fileStream=io.BytesIO(file)
		while True:
			name = fileStream.readline().strip()
			read = fileStream.readline().strip()
			fileStream.readline()
			qual = fileStream.readline().strip()
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
						#print("identical")
						pass
					elif read in s_read:
						#print("read is longer")
						pass
					elif s_read in read:
						#print("s_read is longer")
						pass
					else:
						print("error")
						print(zlib.decompress(zname))
						print(s_read)
						print(read)
				else:
					compressed_data[zname]=(zlib.compress(read),zlib.compress(qual))
					#print(name,read.strip(),qual.strip())
				# check if line is not empty

data={}
				
def gzip_nochunks_byte_u():
	for filename in filenames:
		print("Processing file: ", filename)
		with gzip.open(filename, 'rb') as fh:
			file=fh.read()
			fileStream=io.BytesIO(file)
			while True:
				name = fileStream.readline().strip()
				read = fileStream.readline().strip()
				fileStream.readline()
				qual = fileStream.readline().strip()
				if not (name or read or qual):
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
			fh.close()
		gc.collect()
	#Write the output to disk
	with gzip.open('all.fastq.gz', 'wb') as fh:
		for key, value in data.items():
			fh.write(key+b"\n")
			fh.write(value[0]+b"\n")
			fh.write(b"+")
			fh.write(value[1]+b"\n")
		fh.close()
	#data={}
	gc.collect()
	#Done

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

#ncb=timeit.Timer(gzip_nochunks_byte).timeit(number=1)
ncbu=timeit.Timer(gzip_nochunks_byte_u).timeit(number=1)
#ch=timeit.Timer(gzip_chunks).timeit(number=1)

#print("ncb=",ncb)
print("ncbu=",ncbu)
#print("ch=",ch)

print (time.strftime("%c"))
#input()
