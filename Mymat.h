

#ifndef __MYMAT_H
#define __MYMAT_H

#include <fftw3.h>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string.h>

class Mymat
{

	public:
		Mymat();
		Mymat(int l,int m,int n);
		Mymat(int l,int m,int n,int num);
		Mymat(Mymat &mat1);
		~Mymat();
		
		void rank(int myid, int _size);
		void inposition();
		void outposition();
		void createtype(int n);
		void typefree();
		void createfactor(void);
		void dividefactor();
		void multipfactor();
		void multipfactorx(int k);
		void getplan();
		void fft(int k);
		void ifft(int k);
		void trdFFT(Mymat &mat1, int i, int j, int k);
		void trdIFFT(Mymat &mat1, int i, int j, int k);
		void Laplace(Mymat &mat1, int i, int j, int k);
		void InverseLaplace(Mymat &mat1, int i, int j, int k);
		void NablaTimes(Mymat &V, Mymat &W, Mymat &mat1, int i, int j);
		void Times(Mymat &V, Mymat &W, Mymat &Omega1													, Mymat &Omega2, Mymat &Omega3);	
		void getF(int N);
		void getF1(int N);
		void getU0(int N);
		void getvalue();

		void trans_x(Mymat &mat1);
		void trans_y(Mymat &mat1);
		void trans_z(Mymat &mat1);
		void retrans_x(Mymat &mat1);
		void retrans_y(Mymat &mat1);
		void retrans_z(Mymat &mat1);

		void getVW(Mymat &mat1, Mymat &mat2);

		
		double norm_inf();

		Mymat& operator+=(const Mymat& mat1);	
		Mymat& operator=(const Mymat& mat1);
		Mymat& operator/=(double alpha);
		Mymat operator-(const Mymat& mat1) const;
		Mymat operator+(const Mymat& mat1) const;
		Mymat operator*(double alpha) const;


		double* ele;
		fftw_plan p_dct,p_idct,p_dst,p_idst;
		
		std::vector<double> factor;
	
		void myprint(int a, int b);
	protected:
		//the buff to send
		double* inbuff_x(int i);
		double* inbuff_y(int i);
		double* inbuff_z(int i);
		//the buff to receive
		double* outbuff_x(int i);
		double* outbuff_y(int i);
		double* outbuff_z(int i);






		

	
	private:
		int size,myid,status,plan_on;
		int size_l;
		int size_m;
		int size_n;
		int myorder[3];
		std::vector<int> in_position;
		std::vector<int> out_position_x;
		std::vector<int> out_position_y;
		std::vector<int> out_position_z;
		std::vector<int> xorder;
		std::vector<int> yorder;
		std::vector<int> zorder;
		
		MPI::Datatype byte_type;
		MPI::Datatype tensor1_type;
		MPI::Datatype xtensor0_type;
		MPI::Datatype ycolumn0_type;
		MPI::Datatype ymatrix0_type;
		MPI::Datatype ytensor0_type;
		MPI::Datatype zcolumn0_type;
		MPI::Datatype zmatrix0_type;
		MPI::Datatype ztensor0_type;
		MPI::Datatype utensor2_type;
		MPI::Datatype vcolumn2_type;
		MPI::Datatype vmatrix2_type;
		MPI::Datatype vtensor2_type;
		MPI::Datatype wcolumn2_type;
		MPI::Datatype wmatrix2_type;
		MPI::Datatype wtensor2_type;
};

#endif

