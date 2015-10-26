# Compilers
CC 		= g++

# Objects, final program name
OBJS 	= ssir-mcmc.o timer.o rand.o readdata.o
PROGRAM = ssir-mcmc

# C++ compiler flags
CFLAGS 	= -O2 -w

# Source files with header
%.o: %.cpp %.h
	@$(CC) $(CFLAGS) $< -c -o $@

# Source files without header (should only be main)
%.o: %.cpp
	@$(CC) $(CFLAGS) $< -c -o $@


# Master program compile of objects using NVCC (maybe just CC later)
$(PROGRAM): $(OBJS)
	@$(CC) $(CFLAGS) $^ -o $@ 

# Cleanup
.PHONY: clean
clean:
	@rm -f $(PROGRAM) $(OBJS)