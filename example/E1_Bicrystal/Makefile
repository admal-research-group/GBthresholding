IDIR  += -I../../include/KWC_Simulation

CC= g++

CFLAGS += -Ofast
CFLAGS += -std=c++14
CFLAGS += -lm
CFLAGS += -lfftw3_threads
CFLAGS += -lfftw3

program: 
	$(CC) main.cpp -o main $(CFLAGS) $(IDIR)

