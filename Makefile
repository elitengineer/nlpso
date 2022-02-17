# © 2022 Robert Hoffmann <robert.hoffmann@smail.emt.h-brs.de>
# I'll release this under a license once I decided which.
CC=g++
CFLAGS=-g -Wall -lm

all: circlemap.cpp nlpso.cpp nlpso.hpp
	$(CC) $(CFLAGS) -o build/circlemap circlemap.cpp nlpso.cpp

clean:
	rm build/circlemap