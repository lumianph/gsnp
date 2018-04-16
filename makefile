##################################
## compiler configuration
NVCCFLAGS = -O3 -arch=sm_13


#################################

GSNP_OBJ	= call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o gsnp_kernel.o compress_kernel.o main.o
BANDWIDTH_OBJ	= decompression.o compress_kernel.o
CC              = g++
NVCC            = nvcc
INCLUDE         = -I./cpp -I./inc -I./
LIB             = -L/usr/local/lib -lcutil_x86_64 -lcuda -lcudart
CFLAGS          = -O3 -fomit-frame-pointer -ffast-math -funroll-loops -mmmx -msse -msse2 -msse3 -fmessage-length=0


gsnp: $(GSNP_OBJ)
	$(CC) -o gsnp $(GSNP_OBJ) $(CFLAGS) $(INCLUDE) $(LIB)
decompression: $(BANDWIDTH_OBJ)
	$(CC) -o decompression $(BANDWIDTH_OBJ) $(CFLAGS) $(INCLUDE) $(LIB) 
call_genotype.o: ./cpp/call_genotype.cc
	$(CC) -c ./cpp/call_genotype.cc $(CFLAGS) $(INCLUDE)
chromosome.o: ./cpp/chromosome.cc
	$(CC) -c ./cpp/chromosome.cc $(CFLAGS) $(INCLUDE)
matrix.o: ./cpp/matrix.cc
	$(CC) -c ./cpp/matrix.cc $(CFLAGS) $(INCLUDE)
normal_dis.o: ./cpp/normal_dis.cc
	$(CC) -c ./cpp/normal_dis.cc $(CFLAGS) $(INCLUDE)
prior.o: ./cpp/prior.cc
	$(CC) -c ./cpp/prior.cc $(CFLAGS) $(INCLUDE)
rank_sum.o: ./cpp/rank_sum.cc
	$(CC) -c ./cpp/rank_sum.cc $(CFLAGS) $(INCLUDE)
gsnp_kernel.o: ./kernel/gsnp_kernel.cu
	$(NVCC) -c ./kernel/gsnp_kernel.cu $(NVCCFLAGS) $(INCLUDE)
compress_kernel.o: ./kernel/compress_kernel.cu
	$(NVCC) -c ./kernel/compress_kernel.cu $(NVCCFLAGS) $(INCLUDE)
main.o: ./cpp/main.cc
	$(CC) -c ./cpp/main.cc $(CFLAGS) $(INCLUDE)
decompression.o: ./cpp/decompression.cpp
	$(CC) -c ./cpp/decompression.cpp $(CFLAGS) $(INCLUDE)

all: gsnp decompression

clean:
	rm -f *.o
cleanall: clean
	rm -f gsnp decompression
