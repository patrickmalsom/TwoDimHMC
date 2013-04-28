CC=gcc
#order of arguments: Temp, sigma, sigma^2
CFLAGS=-Wall -Warray-bounds -fopenmp -march=native -O2
LFLAGS=-fopenmp -lm -lgsl

all: TwoDim

TwoDim.o: TwoDim.c
	${CC} ${CFLAGS} ${EXTFLAGS} -c TwoDim.c

TwoDim: TwoDim.o 
	${CC} ${LFLAGS} TwoDim.o -o TwoDim
 
clean: 
	rm -f TwoDim.o TwoDim test.txt

