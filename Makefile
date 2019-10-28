CC = g++

# CPPFLAGS = -g -O3 -Wall -std=c++14
# CPPFLAGS = -g -O3 -funroll-loops -Wno-used-function -std=c++14 -fopenmp -DLIBBSC_OPENMP_SUPPORT
CPPFLAGS = -g -funroll-loops -Wno-used-function -std=c++14 -D_7ZIP_ST
# CPPFLAGS = -O3 -Wno-used-function -std=c++14 -fopenmp -DLIBBSC_OPENMP_SUPPORT

# LIBS = -lpthread -ltbb
LIBS = -pthread -lz #-lm #-ltbb

BSCSRCS = src/libbsc/bsc.cpp src/libbsc/libbsc/adler32/adler32.cpp src/libbsc/libbsc/bwt/divsufsort/divsufsort.c src/libbsc/libbsc/bwt/bwt.cpp src/libbsc/libbsc/coder/coder.cpp src/libbsc/libbsc/coder/qlfc/qlfc.cpp src/libbsc/libbsc/coder/qlfc/qlfc_model.cpp src/libbsc/libbsc/filters/detectors.cpp src/libbsc/libbsc/filters/preprocessing.cpp src/libbsc/libbsc/libbsc/libbsc.cpp src/libbsc/libbsc/lzp/lzp.cpp src/libbsc/libbsc/platform/platform.cpp 
LZMASRC = src/liblzma/lzma.cpp src/liblzma/liblzma/LzmaEnc.c src/liblzma/liblzma/LzmaDec.c src/liblzma/liblzma/Precomp.c src/liblzma/liblzma/7zFile.c src/liblzma/liblzma/Alloc.c src/liblzma/liblzma/LzFind.c src/liblzma/liblzma/7zStream.c

# SRCS = read.cpp hash.cpp compress.cpp fqreader.cpp decompress.cpp main.cpp
# SRCS = src/main.c src/bseq.c src/sketch.c src/util.c $(BSCSRCS) #minimizer_idx.c#hash.cpp
SRCS = src/mstcom.c src/decompress.c src/minimizers.c src/reads.c src/bucket.c src/duplicate.c src/collectnext.c src/output.c src/bseq.c src/sketch.c src/util.c $(BSCSRCS) $(LZMASRC) #minimizer_idx.c#hash.cpp

OBJS = $(SRCS: .c = .o)

EXEC = mstcom


$(EXEC) : $(OBJS)
	$(CC) $(CPPFLAGS) $(LIBS) $^ -o $@

DSRCS = src/decompress.c $(BSCSRCS)

DOBJS = $(DSRCS: .c = .o)

decompress: $(DOBJS)
	$(CC) $(CPPFLAGS) $(LIBS) $^ -o $@

peidx: src/peidx.c 
	$(CC) $(CPPFLAGS) src/peidx.c -o peidx

test1: src/test1.c 
	$(CC) $(CPPFLAGS) src/test1.c -o test1

test2: src/test2.c 
	$(CC) $(CPPFLAGS) src/test2.c -o test2
	
%.o : %.c
	$(CC) -c $(CPPFLAGS) $<

clean:
	rm -f *.o $(EXEC) *.out
