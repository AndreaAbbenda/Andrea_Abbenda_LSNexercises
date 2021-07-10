CC = g++
CFLAGS = -Wall -O3 --std=c++11

PGA : PGA.o random.o
	mpicxx -std=c++0x random.o PGA.o -o PGA

PGA.o : PGA.C random.h
	mpicxx -std=c++0x -c PGA.C -o PGA.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

