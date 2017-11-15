CC=g++
CFLAGS=-Wall -Werror -g
LDFLAGS=-lm -g

EXEC=./exec
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.c=.o)


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

glut.o: glut.cpp glut.hpp matrix.h
	$(CC) -o $@ -cpp $< $(CFLAGS)
	
	
main.o: main.cpp glut.hpp matrix.h
	$(CC) -o $@ -cpp $< $(CFLAGS)

matrix.o: matrix.cpp matrix.h
	$(CC) -o $@ -cpp $< $(CFLAGS)


run: all
	$(EXEC)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
