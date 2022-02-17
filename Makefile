CC=g++
CFLAGS=-Wall -lm

all: circlemap.cpp nlpso.cpp nlpso.hpp
	$(CC) $(CFLAGS) -o build/circlemap circlemap.cpp nlpso.cpp

clean:
	rm build/circlemap