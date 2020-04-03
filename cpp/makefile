CC=g++
CLINKS= -larmadillo -fopenmp
#CFLAGS= -g -std=c++0x -march=native
CFLAGS= -std=c++0x -O3 -march=native

NUM_THREADS=64

all : run

Example : Example.cpp code.o
	$(CC) $(CFLAGS) -o Example Example.cpp code.o $(CLINKS)

Phasespace : Phasespace.cpp code.o numThreads
	$(CC) $(CFLAGS) -o Phasespace Phasespace.cpp code.o $(CLINKS)

Experiment : Experiment.cpp code.o numThreads
	$(CC) $(CFLAGS) -o Experiment Experiment.cpp code.o $(CLINKS)

Convergence : Convergence.cpp code.o numThreads
	$(CC) $(CFLAGS) -o Convergence Convergence.cpp code.o $(CLINKS)

ShieldingFunction : ShieldingFunction.o MicroColony.o
	g++ -g -std=c++0x -o ShieldingFunction ShieldingFunction.o MicroColony.o -O3 -fopenmp

ShieldingFunction.o : ShieldingFunction.cpp
	g++ -g -std=c++0x -c ShieldingFunction.cpp -O3 -fopenmp -larmadillo

ModelContributions : ModelContributions.cpp code.o numThreads
	$(CC) $(CFLAGS) -o ModelContributions ModelContributions.cpp code.o $(CLINKS)

code.o : code.cpp code.hpp
	$(CC) $(CFLAGS) -c code.cpp $(CLINKS)

MicroColony.o : MicroColony.cpp MicroColony.hpp
	g++ -g -std=c++0x -c MicroColony.cpp -O3 -fopenmp -larmadillo

numThreads :
	@echo "int NUM_THREADS = $(NUM_THREADS);" > numThreads.hpp

run : Example Phasespace Experiment Convergence ModelContributions

	mkdir -p data
	mkdir -p data/PhageInfectionProbability

<<<<<<< HEAD
	./Convergence
=======
	#./Convergence
>>>>>>> 880145737b9c4327a084362a992158dc4315db64

	./Example

	./Phasespace

	./Experiment

	./ModelContributions

ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e2.txt)","")
	./ShieldingFunction 100 	>> data/PhageInfectionProbability/gamma_1e2.txt
endif
ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e2_5.txt)","")
	./ShieldingFunction 316		>> data/PhageInfectionProbability/gamma_1e2_5.txt
endif
ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e3.txt)","")
	./ShieldingFunction 1000	>> data/PhageInfectionProbability/gamma_1e3.txt
endif
ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e3_5.txt)","")
	./ShieldingFunction 3162	>> data/PhageInfectionProbability/gamma_1e3_5.txt
endif
ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e4.txt)","")
	./ShieldingFunction 10000	>> data/PhageInfectionProbability/gamma_1e4.txt
endif
ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e4_5.txt)","")
	./ShieldingFunction 31623	>> data/PhageInfectionProbability/gamma_1e4_5.txt
endif
ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e5.txt)","")
	./ShieldingFunction 100000	>> data/PhageInfectionProbability/gamma_1e5.txt
endif
ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e5_5.txt)","")
	./ShieldingFunction 316228	>> data/PhageInfectionProbability/gamma_1e5_5.txt
endif
ifeq ("$(wildcard *data/PhageInfectionProbability/gamma_1e6.txt)","")
	./ShieldingFunction 1000000	>> data/PhageInfectionProbability/gamma_1e6.txt
endif

clean :
	rm -f *.o
	rm -f Example
	rm -f Phasespace
	rm -f Experiment
	rm -f Convergence
	rm -f ModelContributions
	rm -f ShieldingFunction
	rm -f numThreads.hpp
