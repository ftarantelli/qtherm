CC = g++ -g -Wall #-I /home/users/tarante/armadillo-10.8.2/include/

all: therm.out eqth.out

therm.out: qtherm.cpp include/
	$(CC) qtherm.cpp -otherm.out -larmadillo -llapack -lblas -fopenmp
#	./sim.out

eqth.out: eqtherm.cpp include/
	$(CC) eqtherm.cpp -oeqth.out -larmadillo -llapack -lblas -fopenmp
