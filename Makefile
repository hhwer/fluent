
main: test.out
	mpirun -np 8 ./test.out 
test.out: main.o Mymat.o FMymat.o factor.o value.o
	mpic++ -o test.out main.o FMymat.o Mymat.o factor.o value.o Mymat.h -lfftw3  -Wall -std=c++14

main.o: main.cc
	mpic++ -c main.cc -Wall -std=c++14
Mymat.o: Mymat.cc
	mpic++ -c Mymat.cc -Wall -std=c++14
FMymat.o: FMymat.cc
	mpic++ -c FMymat.cc -Wall -std=c++14
factor.o: factor.cc
	mpic++ -c factor.cc -Wall -std=c++14
value.o: value.cc
	mpic++ -c value.cc -Wall -std=c++14



clean:
	rm *.o *.out *.txt
clean1:
	rm  *.txt
gdb: main.cc Mymat.cc FMymat.cc factor.cc value.cc Mymat.h
	mpic++ -g -o test.out main.cc FMymat.cc Mymat.cc factor.cc value.cc Mymat.h -lfftw3  -Wall -std=c++14
run: 
	mpirun -np 8 ./test.out 

t: htest1.out
	mpirun -np 8 ./htest1.out 4
htest1.out: test1.cc Mymat.cc FMymat.cc factor.cc value.cc Mymat.h
	mpic++ -g -o htest1.out test1.cc  FMymat.cc Mymat.cc factor.cc value.cc Mymat.h -lfftw3  -Wall -std=c++14
