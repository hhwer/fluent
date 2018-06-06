
#include "Mymat.h"



/* --------------------------------------------------------------------------*/
/**
* @brief 做3d_fft
*
* @param mat0 数据存贮矩阵
* @param mat1 用来计算FFT的容器
*/
/* ----------------------------------------------------------------------------*/
void Mymat::trdFFT(Mymat &mat1,int i, int j, int k)
{
		this->trans_x(mat1);
		mat1.fft(i);
		this->retrans_x(mat1);
		this->trans_y(mat1);
		mat1.fft(j);
		this->retrans_y(mat1);
		this->trans_z(mat1);
		mat1.fft(k);
		this->retrans_z(mat1);

}


/* --------------------------------------------------------------------------*/
/**
* @brief 做3d_ifft   包含了/N^3
*
* @param mat0 数据存贮矩阵
* @param mat1 用来计算FFT的容器
* @param i,j,k 三个方向上的对称性
*/
/* ----------------------------------------------------------------------------*/
void Mymat::trdIFFT(Mymat &mat1, int i, int j, int k)
{
		this->trans_x(mat1);
		mat1.ifft(i);
		this->retrans_x(mat1);
		this->trans_y(mat1);
		mat1.ifft(j);
		this->retrans_y(mat1);
		this->trans_z(mat1);
		mat1.ifft(k);
		this->retrans_z(mat1);
//除N^3
		double alpha = pow(mat1.size_l*2,3);
		(*this)/=alpha;
//		double* p = &ele[0];
//		for(int i=0;i<size_l*size_m*size_n-1;i++)
//		{
//			*p /= alpha;
//			p++;
//		}
//		*p /= alpha;
}

void Mymat::Laplace(Mymat &mat1, int i, int j, int k)
{
	this->trdFFT(mat1,i,j,k);
	this->multipfactor();
	this->trdIFFT(mat1,i,j,k);
}

void Mymat::InverseLaplace(Mymat &mat1, int i, int j, int k)
{
	this->trdFFT(mat1,i,j,k);
	this->dividefactor();
	this->trdIFFT(mat1,i,j,k);
}


/* --------------------------------------------------------------------------*/
/**
* @brief nable \times u
*
* @param V  用来存储第二个分量
* @param W	用来存储第三个分量 
* @param mat1 做fft的容器
* @param i	V在z上的对称性 1 对称 -1反对称
* @param j	W在y上的的对称性
*
* @returns   
*/
/* ----------------------------------------------------------------------------*/
void Mymat::NablaTimes(Mymat &V, Mymat &W, Mymat &mat1,									int i, int j)
{
	this->getVW(V,W);
	W.trans_y(mat1);
	mat1.fft(j);
	mat1.multipfactorx(j);
	mat1.ifft(-j);
	W.retrans_y(mat1);

	V.trans_z(mat1);
	mat1.fft(i);
	mat1.multipfactorx(i);
	mat1.ifft(-i);
	V.retrans_z(mat1);
	
	(*this) = W-V;
}

/* --------------------------------------------------------------------------*/
/**
* @brief 在U的三个方向补一倍的0
*
* @param U
*/
/* ----------------------------------------------------------------------------*/
void Mymat::bigby(Mymat &U)
{	
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
				ele[i+j*size_l+k*size_l*size_m] = 0;
			}
		}
	}

	for(int k=0;k<U.size_n;k++)
	{
		for(int j=0;j<U.size_m;j++)
		{
			for(int i=0;i<U.size_l;i++)
			{
				ele[i+j*size_l+k*size_l*size_m] = U.ele[i+j*U.size_l+k*U.size_l*U.size_m];
			}
		}
	}


}

void Mymat::smallby(Mymat &bigU)
{	
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
				ele[i+j*size_l+k*size_l*size_m] = 													bigU.ele[i+j*bigU.size_l+k*bigU.size_l*bigU.size_m];
			}
		}
	}
}

void Mymat::Times(Mymat &V, Mymat &W, Mymat &Omega1													, Mymat &Omega2, Mymat &Omega3)
{
	int num;
	this->getVW(V,W);
	Omega1.getVW(Omega2,Omega3);
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
				num = i+j*size_l+k*size_l*size_m;
				ele[num] = Omega2.ele[num]*W.ele[num] 											- Omega3.ele[num]*V.ele[num];
			}
		}
	}

}


/* --------------------------------------------------------------------------*/
/**
* @brief times(omega times u)) 考虑了混淆误差
*
* @param V
* @param W
* @param 
* @param Omega2
* @param Omega3
*/
/* ----------------------------------------------------------------------------*/
void Mymat::Times_v2(Mymat &V, Mymat &W, Mymat &Omega1, Mymat &Omega2, Mymat &Omega3,Mymat &mat1, Mymat &bigU, Mymat &bigV, Mymat &bigW, Mymat &bigOmega1, Mymat &bigOmega2, Mymat &bigOmega3, Mymat &bigmat1)
{
	int num;
	this->trdFFT(mat1,-1,1,1); 
	Omega1.trdFFT(mat1,1,-1,-1);
	bigU.bigby(*this);
	bigOmega1.bigby(Omega1);
	bigU.trdIFFT(bigmat1,-1,1,1);
	bigOmega1.trdIFFT(bigmat1,1,-1,-1);

	bigU.getVW(bigV,bigW);
	bigOmega1.getVW(bigOmega2,bigOmega3);
	for(int k=0;k<bigU.size_n;k++)
	{
		for(int j=0;j<bigU.size_m;j++)
		{
			for(int i=0;i<bigU.size_l;i++)
			{
				num = i+j*bigU.size_l+k*bigU.size_l*bigU.size_m;
				bigU.ele[num] = bigOmega2.ele[num]*bigW.ele[num] 											- bigOmega3.ele[num]*bigV.ele[num];
			}
		}
	}
	bigU.trdIFFT(bigmat1,-1,1,1);
	this->smallby(bigU);
}




/* --------------------------------------------------------------------------*/
/**
* @brief 计算partial Omega / partial t  的半离散的右边 用于RK
*
* @param Omega1
* @param Omega2
* @param Omega3
* @param U
* @param V
* @param W
* @param mat1
* @param nu
* @param tau
*
* @returns   
*/
/* ----------------------------------------------------------------------------*/
double Mymat::f(Mymat &Omega2, Mymat &Omega3, Mymat &U, 						Mymat &V, Mymat &W, Mymat &mat1,Mymat &bigOmega1, Mymat &bigOmega2, Mymat &bigOmega3, Mymat &bigU, Mymat &bigV, Mymat &bigW, Mymat &bigmat1, double nu, double tau, int &sig)
{
	double norm = 0;
	Mymat tempOmega(*this);
	//psi = -laplace^{-1} omega 
	U = (*this);
	U.InverseLaplace(mat1,1,-1,-1);
	U = U*(-1);

	//u(W) = nabla \times psi
	U.NablaTimes(V,W,mat1,-1,-1);  //结果在U中
	if(sig==1)
	{
		sig = U.norm_inf();
		norm = fabs(U.ele[sig]);
	}
	//u = omega \times u    
//	U.Times(V,W,*this,Omega2,Omega3);
	U.Times_v2(V,W,*this,Omega2,Omega3,mat1,bigU,bigV,bigW,bigOmega1,bigOmega2,bigOmega3,bigmat1);


	//U =nabla \times (Omega\time U)
	U.NablaTimes(V,W,mat1,-1,-1);  //结果在U中
	//V 存储 laplace omega
	W = U;
	U = tempOmega;
	U.Laplace(mat1,1,-1,-1);
	*this = U*nu - W;
	*this = (*this)*tau;
	return norm;
}


