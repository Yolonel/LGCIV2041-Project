CC       = gcc  
LD       = gcc
CFLAGS   = -O3 -Dgraphic -Wall -g
LFLAGS   = -Wall -O3 -g
LIBS     = -lm 
#
PROG     = exec
LISTEOBJ = \
  main.o fonction.o matrix_t.o
.c.o :
	$(CC) -c  $(CFLAGS) -o $@ $<

all        : $(PROG)
matrix_t.o : matrix_t.c matrix_t.h
fonction.o : fonction.c fonction.h matrix_t.h
main.o     : main.c fonction.h

$(PROG) : $(LISTEOBJ)
	$(LD) -o $(PROG) $(LFLAGS) $(LISTEOBJ) $(LIBS)

clean :
	rm -vf $(PROG) $(LISTEOBJ) core a.out


