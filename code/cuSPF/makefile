NVCC 	= nvcc
CC 		= g++
NVFLAGS = -arch=sm_20 -O2
CFLAGS 	= -O2
OBJ 	= cuSPF.o timer.o rand.o readdata.o
PROGRAM = cuSPF.x

%.o: %.cpp %.h
	@$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.cu
	@$(NVCC) -c -o $@ $< $(NVFLAGS)

$(PROGRAM): $(OBJ)
	@$(NVCC) -o $@ $^ $(NVFLAGS)

.PHONY: clean

clean:
	@rm -f $(PROGRAM) $(OBJ)