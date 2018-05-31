
main: test
	mpirun -np 8 ./test 64
test: main.o Mymat.o FMymat.o 
	mpic++ -o test main.o FMymat.o Mymat.o Mymat.h -lfftw3  -Wall -std=c++14
main.o: main.cc
	mpic++ -c main.cc -Wall -std=c++14
Mymat.o: Mymat.cc
	mpic++ -c Mymat.cc -Wall -std=c++14
FMymat.o: FMymat.cc
	mpic++ -c FMymat.cc -Wall -std=c++14

clean:
	rm *.o test
gdb: main.cc Mymat.cc FMymat.cc Mymat.h
	mpic++ -o test main.cc FMymat.cc Mymat.cc Mymat.h -lfftw3 -g -Wall -std=c++14
run: 
	mpirun -np 8 ./test 


