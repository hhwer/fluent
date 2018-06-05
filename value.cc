
#include "Mymat.h"


/* --------------------------------------------------------------------------*/
/**
* @brief 生成初值cos x*sin y*sin z
*
* @param N 问题单个方向上的规模
*/
/* ----------------------------------------------------------------------------*/
void Mymat::getF(int N)
{
	for(int k=0;k<size_n;k++)
	{
		double s1 = sin((1.0*(myorder[0]*size_n+k+0.5)/N)*M_PI);
		for(int j=0;j<size_m;j++)
		{
			double s2 = sin((1.0*(myorder[1]*size_m+j+0.5)/N)*M_PI);
			for(int i=0;i<size_l;i++)
			{
				ele[i+j*size_l+k*size_l*size_m] = 												s1*s2*cos((1.0*(myorder[2]*size_l+i+0.5)/N)*M_PI);
			}
		}
	}
}


/* --------------------------------------------------------------------------*/
/**
* @brief 生成初值sin x
*
* @param N 问题单个方向上的规模
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
				ele[i+j*size_l+k*size_l*size_m] = 												sin((1.0*(myorder[2]*size_l+i+0.5)/N)*M_PI);
			}
		}
	}
}
/* --------------------------------------------------------------------------*/
/**
* @brief 生成cos x
*
* @param N 问题单个方向上的规模
*/
/* ----------------------------------------------------------------------------*/
void Mymat::getF2(int N)
{
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
				ele[i+j*size_l+k*size_l*size_m] = 												cos((1.0*(myorder[2]*size_l+i+0.5)/N)*M_PI);
			}
		}
	}
}

/* --------------------------------------------------------------------------*/
/**
* @brief cosy*sinz*sinz*sinx-cosz*siny*siny*sinx
*
* @param N
*/
/* ----------------------------------------------------------------------------*/
void Mymat::getF3(int N)
{
	for(int k=0;k<size_n;k++)
	{
		double s1 = sin((1.0*(myorder[0]*size_n+k+0.5)/N)*M_PI);
		double c1 = cos((1.0*(myorder[0]*size_n+k+0.5)/N)*M_PI);
		for(int j=0;j<size_m;j++)
		{
			double s2 = sin((1.0*(myorder[1]*size_m+j+0.5)/N)*M_PI);
			double c2 = cos((1.0*(myorder[1]*size_m+j+0.5)/N)*M_PI);
			for(int i=0;i<size_l;i++)
			{
				double s3 = sin((1.0*(myorder[2]*size_l+i+0.5)/N)*M_PI);
				ele[i+j*size_l+k*size_l*size_m] = c2*s1*s1*s3-c1*s2*s2*s3;
			}
		}
	}
}
/* --------------------------------------------------------------------------*/
/**
* @brief cosz*sinx*cosy-cos*siny*siny*sinx
*
* @param N
*/
/* ----------------------------------------------------------------------------*/
void Mymat::getF4(int N)
{
	for(int k=0;k<size_n;k++)
	{
		double s1 = sin((1.0*(myorder[0]*size_n+k+0.5)/N)*M_PI);
		double c1 = cos((1.0*(myorder[0]*size_n+k+0.5)/N)*M_PI);
		for(int j=0;j<size_m;j++)
		{
			double s2 = sin((1.0*(myorder[1]*size_m+j+0.5)/N)*M_PI);
			double c2 = cos((1.0*(myorder[1]*size_m+j+0.5)/N)*M_PI);
			for(int i=0;i<size_l;i++)
			{
				double s3 = sin((1.0*(myorder[2]*size_l+i+0.5)/N)*M_PI);
				ele[i+j*size_l+k*size_l*size_m] = c2*s1*s1*s3-c1*s2*s2*s3;
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
				ele[i+j*size_l+k*size_l*size_m] = 												(myorder[2]*size_l+i+1)*(myorder[1]*size_m+j+1)*(myorder[0]*size_n+k+1);
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
			ele[i+j*size_l+k*size_l*size_m] = i*myid;
			}
		}
	}

}



