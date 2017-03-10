# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: ../bin/entropy
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./decimal.f90
main.o: ./main.f90 decimal.o plot.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
plot.o: ./plot.f90 decimal.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plot.f90
util.o: ./util.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
SRC = ./decimal.f90 ./util.f90 ./main.f90 ./plot.f90
OBJ = decimal.o util.o main.o plot.o
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
