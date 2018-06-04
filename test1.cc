
//

#include "Mymat.h"

int main(int argc, char** argv)
{


	MPI::Init(argc, argv);
	int n,Max,myid,totalsize,size;	
	myid = MPI::COMM_WORLD.Get_rank();
	totalsize = MPI::COMM_WORLD.Get_size();
	double tau=1, nu=1;
	int N=pow(2,2);
	if(argc>1)
	{
		N = atoi(argv[1]);
	}

	int aaa =10;
	while(aaa==1);
	{
	}
	size = pow(totalsize+1.0, 1.0/3);
	n = N/size;
	Max = 1;
//	double n3 = pow(N,3);


	Mymat U(n,n,n);
	Mymat V(n,n,n);
	Mymat W(n,n,n);
	Mymat Omega1(n,n,n,0);
	U.rank(myid,size);
	V.rank(myid,size);
	W.rank(myid,size);
	Omega1.rank(myid,size);
	U.createtype(n);
	V.createtype(n);
	W.createtype(n);
	Omega1.createtype(n);
	U.outposition();
	V.outposition();
	W.outposition();
	Omega1.outposition();
	U.createfactor();
	V.createfactor();
	
	Mymat mat1(N,n,n/size);
	mat1.rank(myid,size);
	mat1.inposition();
	mat1.getplan();
	Omega1.getF1(N);
	Omega1.myprint(0,0);

	W=Omega1;
	W.trans_x(mat1);
	mat1.fft(1);
	mat1.multipfactorx(1);
	mat1.ifft(-1);
	W.retrans_x(mat1);
	W.myprint(0,1);
	V.getF2(N);
	V.myprint(0,2);

	U.typefree();
	V.typefree();
	W.typefree();
	Omega1.typefree();

}

