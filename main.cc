//

#include "Mymat.h"

int main(int argc, char** argv)
{
	int n,Max,myid,totalsize,size;
	time_t start,stop,start1,stop1;
	double err[2],total[2];
	double mu=1.5;
	int N=pow(2,4);
	if(argc>1)
	{
		N = atoi(argv[1]);
		if(argc>2)
			mu = atof(argv[2]);
	}



	MPI::Init(argc, argv);
	myid = MPI::COMM_WORLD.Get_rank();
	totalsize = MPI::COMM_WORLD.Get_size();
	size = pow(totalsize+1.0, 1.0/3);
	n = N/size;
	Max = 1000;
	double n3 = pow(N,3);
	start = clock();
	start1 = time(NULL);

	int aaa =0;
	while(aaa==1);
	{
	}


	Mymat mat0(n,n,n,0);
	Mymat mat1(N,n,n/size);
	Mymat U(n,n,n);
	Mymat F(n,n,n);
	mat0.rank(myid,size);
	mat1.rank(myid,size);
	F.rank(myid,size);
	U.rank(myid,size);
	mat0.createtype(n);
	F.getF(N);			//题目的右端项
//	F.getF1(N);         //真解为sinx+icosx时对应的右端项
	mat0.outposition();
	mat1.inposition();
	mat0.createfactor(N,mu);
	
	for(int j=0;j<Max;j++)
    {	
		U = mat0;
		
		mat0^=3;
		mat0 = ((mat0-F)-(U*mu));

		//3d-fft
		mat0.trans_x(mat1);
		mat1.fft();
		mat0.retrans_x(mat1);
		mat0.trans_y(mat1);
		mat1.fft();
		mat0.retrans_y(mat1);
		mat0.trans_z(mat1);
		mat1.fft();

		mat0.retrans_z(mat1);
		mat0.dividefactor();
//		mat0.multipfactor();

		mat0.trans_x(mat1);
		//ifft
		mat1.ifft();

		mat0.retrans_x(mat1);
		mat0.trans_y(mat1);
		//ifft
		mat1.ifft();

		mat0.retrans_y(mat1);
		mat0.trans_z(mat1);
		//ifft
		mat1.ifft();

		mat0.retrans_z(mat1);
		mat0/=n3;
		

		

		



		err[0] = (U-mat0).norm_inf();
		err[1] = U.norm_inf();
		MPI::COMM_WORLD.Allreduce(err, total, 2, MPI_DOUBLE, MPI_MAX);
		if(myid==0)
		{
			std::cout << "relative err= " << total[0]/(total[1]) << " ,err="					<< total[0]	<< " inf_U= "<< total[1] << std::endl;
		}
		
		if(total[0]/total[1]<1e-24 || total[0]==0) 
		{
			stop = clock();
			stop1 = time(NULL);
			if(myid==0)
			{
				std::cout << "cputime:" << 														double(stop-start)/CLOCKS_PER_SEC<< std::endl;
				std::cout << "time:" << 														(stop1-start1)<< std::endl;
				std::cout << "ite:" << j+1 << std::endl;
			}
			break;
		}
  	
	}
	
	mat0.typefree();
	MPI::Finalize();
}

