CC = g++ -g -Wall #-I /home/users/tarante/armadillo-10.8.2/include/

all: therm.out

therm.out: qtherm.cpp include/
	$(CC) qtherm.cpp -otherm.out -larmadillo -llapack -lblas -fopenmp
#	./sim.out
