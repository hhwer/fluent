
#define NN  2760

/**
* @file test1.cc
* @brief 测试函数
* @author hh
* @version 1
* @date 2018-06-06
*/
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

	int aaa =1;
	if(myid==0)
	{
		aaa=0;
	}
//	std::cout<<myid<<aaa<<std::endl;
//	while(aaa==0)
//	{
//	}
	size = pow(totalsize+1.0, 1.0/3);
	n = N/size;
	Max = 1;


	Mymat U(n,n,n);
	Mymat V(n,n,n);
	Mymat W(n,n,n);
	U.rank(myid,size);
	V.rank(myid,size);
	W.rank(myid,size);
	U.createtype(n);
	V.createtype(n);
	W.createtype(n);
	U.outposition();
	V.outposition();
	W.outposition();
	U.createfactor();
	V.createfactor();
	Mymat Omega1(n,n,n,0);
	Omega1.rank(myid,size);
	Omega1.createtype(n);
	Omega1.outposition();
	Mymat Omega2(n,n,n,0);
	Omega2.rank(myid,size);
	Omega2.createtype(n);
	Omega2.outposition();
	Mymat Omega3(n,n,n,0);
	Omega3.rank(myid,size);
	Omega3.createtype(n);
	Omega3.outposition();
	
	Mymat mat1(N,n,n/size);
	mat1.rank(myid,size);
	mat1.inposition();
	mat1.getplan();

if(myid==0){
	std::ifstream infile;
	std::ofstream outfile;
	infile.open("normW.txt",std::ios::in);
	
	if(!infile.is_open())
		std::cout << "open infile failure" << std::endl;

	double num[NN];
	int i =0;
	while(!infile.eof())
	{
		infile >> num[i];
		i++;
	}
	infile.close();
	
	outfile.open("logW.txt",std::ios::app);
	if(!outfile.is_open())
		std::cout << "open file failure" << std::endl;
	for(int i =0;i<NN-1;i++)
	{
		outfile << log(num[i+1])-log(num[i]) << std::endl;		
	}
	outfile.close();
	
	outfile.open("slopeW.txt",std::ios::app);
	if(!outfile.is_open())
		std::cout << "open file failure" << std::endl;
	for(int i =0;i<NN-1;i++)
	{
		outfile << num[i+1]-num[i] << std::endl;		
	}
	outfile.close();
	infile.open("normU.txt",std::ios::in);
	
	if(!infile.is_open())
		std::cout << "open infile failure" << std::endl;

	i=0;
	while(!infile.eof())
	{
		infile >> num[i];
		i++;
	}
	infile.close();
	
	outfile.open("logU.txt",std::ios::app);
	if(!outfile.is_open())
		std::cout << "open file failure" << std::endl;
	for(int i =0;i<NN-1;i++)
	{
		outfile << log(num[i+1])-log(num[i]) << std::endl;		
	}
	outfile.close();
	
	outfile.open("slopeU.txt",std::ios::app);
	if(!outfile.is_open())
		std::cout << "open file failure" << std::endl;
	for(int i =0;i<NN-1;i++)
	{
		outfile << num[i+1]-num[i] << std::endl;		
	}
	outfile.close();
}
//	U.getF(N);
//	U.myprint(0,0);
//	int num = U.norm_inf();
//	double norm = fabs(U.ele[num]);
//	std::cout << num << ' ' << norm << std::endl;

//// Nabla \times u
//	U.getF(N);
//	U.NablaTimes(V,W,mat1,-1,-1);
//	W.myprint(0,1);


//// omega times u
//	Omega1.getF(N);
//	U.getF1(N);
//	U.Times(V,W,Omega1,Omega2,Omega3);
//	V.getF3(N);
//	U.myprint(0,1);
//	V.myprint(0,2);


//f= cosx*siny*sinz  laplace(f)=-3f
//	U.getF(N);
//	V = U*(-3);
//	U.Laplace(mat1,1,-1,-1);
//	U.myprint(0,1);
//	V.myprint(0,2);

////f= cosx*siny*sinz  inverselaplace(f)=-f/3
//	U.getF(N);
//	V = U*(-1.0/3);
//	V.myprint(0,2);
//	U.InverseLaplace(mat1,1,-1,-1);
//	U.myprint(0,1);


////i=1,-1  u=cosx,sinx 求导
//	int i =2;  
//	if(i==1)
//	{
//		Omega1.getF2(N);
//		V.getF1(N);
//	}
//	else
//	{
//		Omega1.getF1(N);
//		V.getF2(N);
//	}
//	Omega1.myprint(0,0);
//	W=Omega1;
//	W.trans_x(mat1);
//	mat1.fft(i);
//	mat1.multipfactorx(i);
//	mat1.ifft(-i);
//	W.retrans_x(mat1);
//	W.myprint(0,1);
//	V=V*(-i);	
//	W.myprint(0,2);


	U.typefree();
	V.typefree();
	W.typefree();
	Omega1.typefree();
	Omega2.typefree();
	Omega3.typefree();






	MPI::Finalize();

}

