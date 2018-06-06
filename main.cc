//

#include "Mymat.h"

int main(int argc, char** argv)
{


	MPI::Init(argc, argv);
	
	clock_t start,stop,start1,stop1;
	start = clock();
	start1 = time(NULL);

	int n,Max,myid,totalsize,size,num0,num1,num2;	
	double norm1,norm2;
	myid = MPI::COMM_WORLD.Get_rank();
	totalsize = MPI::COMM_WORLD.Get_size();
	double tau=0.01, nu=0;
	int N=pow(2,2);
	Max = 1000;
	if(argc>1)
	{
		N = atoi(argv[1]);
		if(argc>2)
			Max = atoi(argv[2]);
	}

	int aaa =10;
	
	if(myid==0)
		aaa=11;

	while(aaa==1)
	{
	}
	size = pow(totalsize+1.0, 1.0/3);
	n = N/size;
	int *gathernum;
	double *gathernorm;
	std::vector<double> v(totalsize);
	if(myid==0)
	{
		gathernum = (int*) malloc(sizeof(int)*totalsize);
		gathernorm = (double*)malloc(sizeof(double)*totalsize);
		
	}
	std::ofstream fp;


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

	

	Omega1.getOmega0(N);

	for(int ite=0;ite<Max;ite++)
	{	
		num1 = 1;
		K0 = Omega1;
		norm1 = Omega1.f(Omega2, Omega3, U, V, W, mat1, nu, tau,num1);
		K1 = Omega1;
		Omega1 = K0 + K1*0.5;
		Omega1.f(Omega2, Omega3, U, V, W, mat1, nu, tau,N);
		K2 = Omega1;
		Omega1 = K0 + K2*0.5;
		Omega1.f(Omega2, Omega3, U, V, W, mat1, nu, tau,N);
		K3 = Omega1;
		Omega1 += K0;
		Omega1.f(Omega2, Omega3, U, V, W, mat1, nu, tau,N);
		Omega1 = K0 + (K1+K2*2+K3*2+Omega1)*(1.0/6);
		num2 = Omega1.norm_inf();
		norm2 = fabs(Omega1.ele[num2]);	
		if(myid==0)
		std::cout<<ite<<std::endl;

		MPI_Gather(&num1, 1, MPI::INT,gathernum,1,MPI::INT										,0,MPI::COMM_WORLD);
		MPI_Gather(&norm1, 1, MPI::DOUBLE,gathernorm											,1,MPI::DOUBLE,0,MPI::COMM_WORLD);
		if(myid==0)
		{	
			for(int i=0;i<totalsize;i++) 
			{
				v[i] = gathernorm[i];
			}
			num0 = (int) (std::max_element(v.begin(),v.end()) - v.begin());
			num1 = gathernum[num0];
			norm1 = gathernorm[num0];
			fp.open("U.txt",std::ios::app);
			fp.setf(std::ios::fixed);
			fp<<std::setprecision(6)<<(num0%size*n+num1%n+0.5)*M_PI/N;
			fp<<' ';
			fp<<std::setprecision(6)<<(num0/size%size*n+num1/n%n+0.5)*M_PI/N;
			fp<<' ';
			fp<<std::setprecision(6)<<(num0/size/size*n+num1/n/n+0.5)*M_PI/N<<std::endl;
			fp.close();
			fp.open("normU.txt",std::ios::app);
			fp<<norm1<<std::endl;
//			fp<<log(norm1)<<std::endl;
			fp.close();
		}
		
//		MPI_Gather(&num2, 1, MPI::INT,gathernum,1,MPI::INT										,0,MPI::COMM_WORLD);
//		MPI_Gather(&norm2, 1, MPI::DOUBLE,gathernorm											,1,MPI::DOUBLE,0,MPI::COMM_WORLD);
//		if(myid==0)
//		{	
//			for(int i=0;i<totalsize;i++) 
//			{
//				v[i] = gathernorm[i];
//			}
//			num0 = (int) (std::max_element(v.begin(),v.end()) - v.begin());
//			num1 = gathernum[num0];
//			norm1 = gathernorm[num0];
//			fp.open("W.txt",std::ios::app);
//			fp.setf(std::ios::fixed);
//			fp<<std::setprecision(6)<<(num0%size*n+num1%n+0.5)*M_PI/N;
//			fp<<' ';
//			fp<<std::setprecision(6)<<(num0/size%size*n+num1/n%n+0.5)*M_PI/N;
//			fp<<' ';
//			fp<<std::setprecision(6)<<(num0/size/size*n+num1/n/n+0.5)*M_PI/N<<std::endl;
//			fp.close();
//			fp.open("normW.txt",std::ios::app);
//			fp<<norm1<<std::endl;
////			fp<<log(norm1)<<std::endl;
//			fp.close();
//		}

	}
//		
//
	stop = clock();
	stop1 = time(NULL);
	if(myid==0)
	{
		std::cout << "cputime:" << 														double(stop-start)/CLOCKS_PER_SEC<< std::endl;
		std::cout << "time:" << 														(stop1-start1)<< std::endl;
	}
	U.typefree();
	V.typefree();
	W.typefree();
	Omega1.typefree();
	Omega2.typefree();
	Omega3.typefree();
	std::vector<double>().swap(v);

	MPI::Finalize();

}

