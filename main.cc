//

#include "Mymat.h"

int main(int argc, char** argv)
{


	MPI::Init(argc, argv);
	
	clock_t start,stop,start1,stop1;
	start = clock();
	start1 = time(NULL);

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
	while(aaa==1)
	{
	}
	size = pow(totalsize+1.0, 1.0/3);
	n = N/size;
	Max = 10;


	Mymat U(n,n,n);
	Mymat V(n,n,n);
	Mymat W(n,n,n);
	Mymat K0(n,n,n);
	Mymat K1(n,n,n);
	Mymat K2(n,n,n);
	Mymat K3(n,n,n);
	Mymat Omega1(n,n,n,0);
	Mymat Omega2(n,n,n);
	Mymat Omega3(n,n,n);
	U.rank(myid,size);
	V.rank(myid,size);
	W.rank(myid,size);
	Omega1.rank(myid,size);
	Omega2.rank(myid,size);
	Omega3.rank(myid,size);
	U.createtype(n);
	V.createtype(n);
	W.createtype(n);
	Omega1.createtype(n);
	Omega2.createtype(n);
	Omega3.createtype(n);
	U.outposition();
	V.outposition();
	W.outposition();
	Omega1.outposition();
	Omega2.outposition();
	Omega3.outposition();
	U.createfactor();
	V.createfactor();
	
	Mymat mat1(N,n,n/size);
	mat1.rank(myid,size);
	mat1.inposition();
	mat1.getplan();

	

	Omega1.getF(N);

	for(int j=0;j<Max;j++)
	{	
		K0 = Omega1;
		K1 = f(Omega1, Omega2, Omega3, U, V, W, mat1, nu, tau);
		Omega1 = K0 + K1*0.5;
		K2 = f(Omega1, Omega2, Omega3, U, V, W, mat1, nu, tau);
		Omega1 = K0 + K2*0.5;
		K3 = f(Omega1, Omega2, Omega3, U, V, W, mat1, nu, tau);
		Omega1 = K0 + K3;
		Omega1 = f(Omega1, Omega2, Omega3, U, V, W, mat1, nu, tau);
		Omega1 = K0 + (K1+K2*2+K3*2+Omega1)*(1.0/6);
		if(myid==0)
		std::cout<<j<<std::endl;
	}
//		
//
//		
//
//		
//
//
//
//		err[0] = (U-mat0).norm_inf();
//		err[1] = U.norm_inf();
//		MPI::COMM_WORLD.Allreduce(err, total, 2, MPI_DOUBLE, MPI_MAX);
//		if(myid==0)
//		{
//			std::cout << "relative err= " << total[0]/(total[1]) << " ,err="					<< total[0]	<< " inf_U= "<< total[1] << std::endl;
//		}/
//		
//		if(total[0]/total[1]<1e-24 || total[0]==0) 
//		{
			stop = clock();
			stop1 = time(NULL);
			if(myid==0)
			{
				std::cout << "cputime:" << 														double(stop-start)/CLOCKS_PER_SEC<< std::endl;
				std::cout << "time:" << 														(stop1-start1)<< std::endl;
			}
//			break;
//		}
//  	
//	}
	U.typefree();
	V.typefree();
	W.typefree();
	Omega1.typefree();
	Omega2.typefree();
	Omega3.typefree();


	MPI::Finalize();

}

