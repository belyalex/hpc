CC = g++
ccflags = -c -std=c++17 -g

all: MulMatr

MulMatr: main.o matrix.o
	$(CC) main.o matrix.o -o MulMatr

main.o: main.cpp matrix.h timer.h
	$(CC) $(ccflags) main.cpp

matrix.o: matrix.cpp matrix.h
	$(CC) $(ccflags) matrix.cpp

clean:
	rm -f *.o MulMatr
