# Compilers
CC 		= g++

# Objects, final program name
OBJS 	= ssir.o timer.o rand.o
PROGRAM = ssir

# C++ compiler flags
CFLAGS 	= -O2 -w

# Source files with header
%.o: %.cpp %.h
	@$(CC) $(CFLAGS) $< -c -o $@

# Source files without header (should only be main.cpp)
%.o: %.cpp
	@$(CC) $(CFLAGS) $< -c -o $@


# Master program compile of objects using NVCC (maybe just CC later)
$(PROGRAM): $(OBJS)
	@$(CC) $(CFLAGS) $^ -o $@ 

# Cleanup
.PHONY: clean
clean:
	@rm -f $(PROGRAM) $(OBJS)