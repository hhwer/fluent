
/**
* @file factor.cc
* @brief 
* @author hh
* @version hh
* @date 2018-06-01
*/

#include "Mymat.h"

/* --------------------------------------------------------------------------*/
/**
* @brief 生成每个数据块中的laplace对应的傅里叶系数-i^2-j^2-k^2
*/
/* ----------------------------------------------------------------------------*/
void Mymat::createfactor(void)
{	
	factor.resize(size_l*size_m*size_n);
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
				factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i),2);
			}
		}
	}

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
	}
	*p *= *f;
	p++;
	*p *= *f;
}

/* --------------------------------------------------------------------------*/
/**
* @brief 从u的傅里叶系数变成partial_y u的系数
*/
/* ----------------------------------------------------------------------------*/
void Mymat::multipfactory(void)
{	
	int num = 0;
	double a;
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			int alpha = myorder[1]*size_m+j;
			for(int i=0;i<size_l;i++)
			{	
				a = alpha*ele[num][0];
				ele[num][0] = -alpha*ele[num][1];
				ele[num][1] = a;
				num += 1;
			}
		}
	}
}

/* --------------------------------------------------------------------------*/
/**
* @brief 从u的傅里叶系数变成partial_z u的系数
*/
/* ----------------------------------------------------------------------------*/
void Mymat::multipfactorz(void)
{	
	int num = 0;
	double a;
	for(int k=0;k<size_n;k++)
	{
		int alpha = myorder[0]*size_n+k;
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{	
				a = alpha*ele[num][0];
				ele[num][0] = -alpha*ele[num][1];
				ele[num][1] = a;
				num += 1;
			}
		}
	}
}

