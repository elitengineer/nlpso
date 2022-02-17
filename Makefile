CC=g++
CFLAGS=-Wall -lm

all: circlemap.cpp nlpso.cpp nlpso.hpp
	$(CC) $(CFLAGS) -o circlemap circlemap.cpp nlpso.cpp

clean:
	rm circlemap