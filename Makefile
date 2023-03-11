CC = g++ -g -Wall #-I /home/users/tarante/armadillo-10.8.2/include/

all: therm.out

therm.out: qtherm.cpp qsystem.hpp
	$(CC) qtherm.cpp -otherm.out -larmadillo -llapack -lblas
#	./sim.out
