
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
		double alpha = pow(mat1.size_l*2-2,3);
		double* p = &ele[0][0];
		for(int i=0;i<2*size_l*size_m*size_n-1;i++)
		{
			*p /= alpha;
			p++;
		}
		*p /= alpha;
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
* @param W	用来存储第三个分量 并且存储结果
* @param mat1 做fft的容器
* @param i	V在z上的对称性 1 对称 -1反对称
* @param j	W在y上的的对称性
*
* @returns   
*/
/* ----------------------------------------------------------------------------*/
void Mymat::NablaTimes(Mymat &V, Mymat &W, Mymat &mat1,									int i, int j)
{
	double alpha = mat1.size_l*2-2;
	this->getVW(V,W);
	W.trans_y(mat1);
	mat1.fft(j);
	W.retrans_y(mat1);
	W.multipfactory();
	W.trans_y(mat1);
	mat1.ifft(-j);
	W.retrans_y(mat1);

	V.trans_z(mat1);
	mat1.fft(i);
	V.retrans_z(mat1);
	V.multipfactorz();
	V.trans_z(mat1);
	mat1.ifft(-i);
	V.retrans_z(mat1);
	
	W= W-V;
	double* p = &W.ele[0][0];
	for(int i=0;i<2*size_l*size_m*size_n-1;i++)
	{
		*p /= alpha;
		p++;
	}
	*p /= alpha;
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
				ele[num][0] = Omega2.ele[num][0]*W.ele[num][0] 					- Omega2.ele[num][1]*W.ele[num][1] - Omega3.ele[num][0]*V.ele[num][0] 		+ Omega3.ele[num][1]*V.ele[num][1];
				ele[num][1] = Omega2.ele[num][0]*W.ele[num][1] 					+ Omega2.ele[num][1]*W.ele[num][0] - Omega3.ele[num][0]*V.ele[num][1] 		- Omega3.ele[num][1]*V.ele[num][0];
			}
		}
	}




}

/* --------------------------------------------------------------------------*/
/**
* @brief 生成右端项
*
* @param N 问题单个方向上的规模
*/
/* ----------------------------------------------------------------------------*/
void Mymat::getF(int N)
{
	for(int k=0;k<size_n;k++)
	{
		double c1 = cos((2.0*(myorder[0]*size_n+k)/N-1)*M_PI);
		double s1 = sin((2.0*(myorder[0]*size_n+k)/N-1)*M_PI);
		for(int j=0;j<size_m;j++)
		{
			double c2 = cos((2.0*(myorder[1]*size_m+j)/N-1)*M_PI);
			double s2 = sin((2.0*(myorder[1]*size_m+j)/N-1)*M_PI);
			for(int i=0;i<size_l;i++)
			{
				ele[i+j*size_l+k*size_l*size_m][0] = 												s1*s2*sin((2.0*(myorder[2]*size_l+i)/N-1)*M_PI);
				ele[i+j*size_l+k*size_l*size_m][1] = 												c1*c2*cos((2.0*(myorder[2]*size_l+i)/N-1)*M_PI);
			}
		}
	}

}

/* --------------------------------------------------------------------------*/
/**
* @brief F=sin(x)+icos(x)时的右端项
*
* @param N
*/
/* ----------------------------------------------------------------------------*/
void Mymat::getF1(int N)
{
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
				ele[i+j*size_l+k*size_l*size_m][0] = 												2 * sin((2*(myorder[2]*size_l+i+0.0)/N-1)*M_PI);
				ele[i+j*size_l+k*size_l*size_m][1] = 												2 * cos((2*(myorder[2]*size_l+i+0.0)/N-1)*M_PI);
			}
		}
	}
}


void Mymat::getU0(int N)
{
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
				ele[i+j*size_l+k*size_l*size_m][0] = 												(myorder[2]*size_l+i+1)*(myorder[1]*size_m+j+1)*(myorder[0]*size_n+k+1);
			}
		}
	}
}
void Mymat::getvalue(void)
{
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
			ele[i+j*size_l+k*size_l*size_m][0] = i*myid;
			}
		}
	}

}


