# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: ../bin/entropy
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./decimal.f90
example_burger.o: ./example_burger.f90 decimal.o plot.o tipos.o numeric_schemes.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./example_burger.f90
main.o: ./main.f90 example_burger.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
numeric_schemes.o: ./numeric_schemes.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./numeric_schemes.f90
plot.o: ./plot.f90 decimal.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plot.f90
tipos.o: ./tipos.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./tipos.f90
util.o: ./util.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
SRC = ./main.f90 ./example_burger.f90 ./util.f90 ./decimal.f90 ./numeric_schemes.f90 ./tipos.f90 ./plot.f90
OBJ = main.o example_burger.o util.o decimal.o numeric_schemes.o tipos.o plot.o
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
