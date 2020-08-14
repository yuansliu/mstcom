# mstcom

The mstcom program is written in C++11 and works on Linux. It is availble under an open-source license.

## Download & Install

	git clone https://github.com/yuansliu/mstcom.git
 	cd mstcom
	make
  
 ## Usage
 To compress:
 
 	./mstcom e -i IN.fastq -o OUTPUT
 	./mstcom e -i IN.fastq -o OUTPUT -p                     #order-preserving mode
 	./mstcom e -i IN_1.fastq -f IN_2.fastq -o OUTPUT
 	./mstcom e -i IN_1.fastq -f IN_2.fastq -o OUTPUT -p     #order-preserving mode
  
 To decompress:
 
 	./mstcom d -i Compressed_Result -o OUTPUT
  
Options

 	-h              print help message
 	-t              number of threads, default: 24
  
