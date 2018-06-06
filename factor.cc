
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
* @brief 生成每个数据块中的laplace对应的傅里叶系数-i^2-(j+1)^2-(k+1)^2
*	仅对应x对称,yz反对称的情形
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
				factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j+1),2);
			}
		}
	}

	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
				factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k+1),2);
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
	double* p = &ele[0];
	double* f = &factor[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
	{
		*p /= *f;
		p++;
		f++;
	}
	*p /= *f;
}


void Mymat::multipfactor(void)
{
	double* p = &ele[0];
	double* f = &factor[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
	{
		*p *= *f;
		p++;
		f++;
	}
	*p *= *f;
}


/* --------------------------------------------------------------------------*/
/**
* @brief 从u的傅里叶系数变成partial_x u的系数 并完成位移 
*以及归一化操作 /N
*/
/* ----------------------------------------------------------------------------*/
void Mymat::multipfactorx(int k)
{	
	int num = -1;
	double alpha = 2.0*size_l;
	if(k==1)
	{
    	for(int j=0;j<size_m*size_n;j++)
	    {	
			num++;
	    	for(int i=1;i<size_l;i++)
	    	{	
				//ele[num] = -i*ele[num++];
	    		
				ele[num] = -i*ele[num+1]/alpha;
				num++;
			}
			ele[num] = 0;
		}
	}
	else
	{	
		num = size_m*size_n*size_l;
    	for(int j=0;j<size_m*size_n;j++)
	    {
			num--;
	    	for(int i=size_l-1;i>0;i--)
	    	{	
	    		ele[num] = i*ele[num-1]/alpha;
				num--;
			}
			ele[num] = 0;
		}
	}
}

