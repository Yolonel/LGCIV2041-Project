CC=g++
CFLAGS=-Wall -Werror -g
LDFLAGS=-lm -g

EXEC=./exec
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.c=.o)


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

fonction.o: fonction.cpp fonction.hpp matrix_t.h
	$(CC) -o $@ -cpp $< $(CFLAGS)
	
	
main.o: main.cpp fonction.hpp matrix_t.h
	$(CC) -o $@ -cpp $< $(CFLAGS)

matrix_t.o: matrix_t.cpp matrix_t.h
	$(CC) -o $@ -cpp $< $(CFLAGS)


run: all
	$(EXEC)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
