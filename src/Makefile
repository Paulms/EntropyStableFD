# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: ../bin/entropy
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./decimal.f90
example_burger.o: ./example_burger.f90 decimal.o plot.o tipos.o numeric_schemes.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./example_burger.f90
main.o: ./main.f90 example_burger.o test2.o test3.o test4.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
numeric_schemes.o: ./numeric_schemes.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./numeric_schemes.f90
numeric_schemes_nd.o: ./numeric_schemes_nd.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./numeric_schemes_nd.f90
plot.o: ./plot.f90 decimal.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plot.f90
test2.o: ./test2.f90 decimal.o plot.o tipos.o numeric_schemes.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test2.f90
test3.o: ./test3.f90 decimal.o plot.o tipos.o numeric_schemes.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test3.f90
test4.o: ./test4.f90 decimal.o plot.o tipos.o numeric_schemes_nd.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test4.f90
tipos.o: ./tipos.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./tipos.f90
util.o: ./util.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
SRC = ./numeric_schemes.f90 ./main.f90 ./test3.f90 ./numeric_schemes_nd.f90 ./test2.f90 ./decimal.f90 ./test4.f90 ./tipos.f90 ./example_burger.f90 ./plot.f90 ./util.f90
OBJ = numeric_schemes.o main.o test3.o numeric_schemes_nd.o test2.o decimal.o test4.o tipos.o example_burger.o plot.o util.o
clean: neat
	-rm -f .../bin/entropy.cppdefs $(OBJ) ../bin/entropy
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
../bin/entropy: $(OBJ) 
	$(LD) $(OBJ) -o ../bin/entropy  $(LDFLAGS)
