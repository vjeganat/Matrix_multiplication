# Before make, you must include
# $ module load gcc/4.9.3 mkl

# If you want to use Intel compilers, you need to modify the following
# flags to link with MKL. 

CC = gcc 
OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT) -m64 -I${MKLROOT}/include -mavx -ftree-vectorize  -ffast-math -fopt-info-vec-optimized  -fpeel-loops #-funsafe-math-optimizations
LDFLAGS = -Wall
# librt is needed for clock_gettime
LDLIBS = -lrt    -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

targets = benchmark-naive benchmark-blocked benchmark-blas
objects = benchmark.o dgemm-naive.o dgemm-blocked.o dgemm-blas.o  

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

benchmark-naive : benchmark.o dgemm-naive.o 
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blocked : benchmark.o dgemm-blocked.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blas : benchmark.o dgemm-blas.o
	$(CC) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects)
