
#include "Mymat.h"


/* --------------------------------------------------------------------------*/
/**
* @brief 从\Delta U-\mu U的fourier系数得到U的傅里叶系数所需的除数
*
* @param N  单个方向的规模
* @param mu	
*/
/* ----------------------------------------------------------------------------*/
void Mymat::createfactor(int N,double mu)
{
	factor.resize(size_l*size_m*size_n);
	if(2*myorder[2]<size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i),2) - mu;
				}
			}
		}
	}
	else if(2*myorder[2]>size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i-N),2) - mu;
				}
			}
		}

	}
	else
	{	
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l/2;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i),2) - mu;
				}
				for(int i=size_l/2;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i-N),2) - mu;
				}
			}
		}
	 }
	

	if(2*myorder[1]<size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j),2);
				}
			}
		}
	}
	else if(2*myorder[1]>size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j-N),2);
				}
			}
		}

	}
	else
	{	
		for(int k=0;k<size_n;k++)
		{
			for(int i=0;i<size_l;i++)
			{
				for(int j=0;j<size_m/2;j++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j),2);
				}
				for(int j=size_m/2;j<size_m;j++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j-N),2);
				}
			}
		}
	}

	if(2*myorder[0]<size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k),2);
				}
			}
		}
	}
	else if(2*myorder[0]>size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k-N),2);
				}
			}
		}

	}
	else
	{	
		for(int i=0;i<size_l;i++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int k=0;k<size_n/2;k++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k),2);
				}
				for(int k=size_n/2;k<size_n;k++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k-N),2);
				}
			}
		}
	}
}

/* --------------------------------------------------------------------------*/
/**
* @brief 从\Delta U -\mu U的fourier得到U的fourier
*/
/* ----------------------------------------------------------------------------*/
void Mymat::dividefactor(void)
{
	if (factor[0]==0)
	{
		ele[0][0] = 0.0;
		ele[0][1] = 0.0;
	}
	else
	{
		ele[0][0] /= factor[0];
		ele[0][1] /= factor[0];
	}
	double* p = &ele[1][0];
	double* f = &factor[1];
	for(int i=1;i<size_l*size_m*size_n-1;i++)
	{
		*p /= *f;
		p++;
		*p /= *f;
		p++;
		f++;
//		ele[i][0] /= factor[i];
//		ele[i][1] /= factor[i];
	}
	*p /= *f;
	p++;
	*p /= *f;
}


void Mymat::multipfactor(void)
{
	double* p = &ele[0][0];
	double* f = &factor[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
	{
		*p *= *f;
		p++;
		*p *= *f;
		p++;
		f++;
//		ele[i][0] /= factor[i];
//		ele[i][1] /= factor[i];
	}
	*p *= *f;
	p++;
	*p *= *f;
}

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
* @brief 做3d_ifft
*
* @param mat0 数据存贮矩阵
* @param mat1 用来计算FFT的容器
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


