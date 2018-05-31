/**
* @file Mymat.cc
* @brief class Mymat implement
* @author hh
* @version 0
* @date 2017-12-09
*/
//


#include "Mymat.h"

Mymat::Mymat(void)
{
	size_l = 0;
	size_m = 0;
	size_n = 0;
	status = 0;
}

Mymat::Mymat(int l, int m, int n )
{
	size_l = l;
	size_m = m;
	size_n = n;
	status = 0;
	ele = (fftw_complex*) malloc( sizeof(fftw_complex)*l*m*n );
	double* p=&ele[0][0];
	for(int i=0;i<2*l*m*n-1;i++)
	{
		*p = 0.0;
		p++;
	}
	*p = 0.0;
}


Mymat::Mymat(int l, int m, int n,int num)
{
	size_l = l;
	size_m = m;
	size_n = n;
	status = 0;
	double num0 = (double)num; 
	ele = (fftw_complex*) malloc( sizeof(fftw_complex)*l*m*n );
	double* p=&ele[0][0];
	for(int i=0;i<2*l*m*n-1;i++)
	{
		*p = num0;
		p++;
	}
	*p = num0;
}


Mymat::Mymat(Mymat& mat1)
{
	size_l = mat1.size_l;
	size_m = mat1.size_m;
	size_n = mat1.size_n;
	status = 0;
	ele = (fftw_complex*) malloc( sizeof(fftw_complex)*size_l*size_m*size_n );
	double* p=&ele[0][0];
	for(int i=0;i<2*size_l*size_m*size_n-1;i++)
	{
		*p = 0.0;
		p++;
	}
}
Mymat::~Mymat(void)
	{
		free(ele);
		if(status)
		{
			typefree();
		}
	}

/* --------------------------------------------------------------------------*/
/**
* @brief 计算数据块在整个方体中所处的位置,并以其为起始,以增序排列各方向上的进程号
*
* @param _myid
* @param _size
*/
/* ----------------------------------------------------------------------------*/
void Mymat::rank(int _myid, int _size)
{
	myid = _myid;
	size = _size;
	xorder.resize(size);
	yorder.resize(size);
	zorder.resize(size);

	myorder[2] = myid%size;
	myorder[0] = myid/size/size;
	myorder[1] = myid/size%size;
	for(int i=0;i<size;i++)
	{
		xorder[i] = (myorder[2]+i)%size + myorder[1] * size +												 myorder[0] * size * size;
		yorder[i] = myorder[2] + (myorder[1]+i)%size * size +												 myorder[0] * size * size;
		zorder[i] = myorder[2] + myorder[1] * size +												 (myorder[0]+i)%size * size * size;
	}
	
}



/* --------------------------------------------------------------------------*/
/**
* @brief 用来计算fft的内存块中,第i块的起始地址
*/
/* ----------------------------------------------------------------------------*/
void Mymat::inposition(void)
{
	in_position.resize(size);
	for(int i=0;i<size;i++)
	{
		in_position[i] = i*size_l/size; 
	}
}

/* --------------------------------------------------------------------------*/
/**
* @brief 计算存储数据的内存块内,三个方向分别剖分成size份时,第i份开始的位置
*/
/* ----------------------------------------------------------------------------*/
void Mymat::outposition(void)
{
	out_position_x.resize(size);
	out_position_y.resize(size);
	out_position_z.resize(size);
	
	for(int i=0;i<size;i++)
	{
		out_position_x[i] = i*size_l*size_m*size_n/size; 
		out_position_y[i] = i*size_l*size_m*size_n/size; 
		out_position_z[i] = i*size_l*size_m/size; 
	}
}

/* --------------------------------------------------------------------------*/
/**
* @brief 生成用于mpi传输的各类datatype
*
* @param n 每个数据块的尺寸
*/
/* ----------------------------------------------------------------------------*/
void Mymat::createtype(int n)
{
	status = 1;
	 byte_type = MPI::DOUBLE.Create_vector(1, 2, 2);
	byte_type.Commit();

	 tensor1_type = byte_type.Create_vector(n*n/size, n, n*size);
	tensor1_type.Commit();

	 xtensor0_type = byte_type.Create_vector(1, n*n*n/size, 0);
	xtensor0_type.Commit();

	 ycolumn0_type = byte_type.Create_vector(n, 1, n);
	ycolumn0_type.Commit();

	 ymatrix0_type = ycolumn0_type.Create_hvector(n, 1, sizeof(fftw_complex));
	ymatrix0_type.Commit();

	 ytensor0_type = ymatrix0_type.Create_vector(n/size, 1, 1);
	ytensor0_type.Commit();
	
	
	 zcolumn0_type = byte_type.Create_vector(n, 1, n*n);
	zcolumn0_type.Commit();

	 zmatrix0_type = zcolumn0_type.Create_hvector(n, 1, 											sizeof(fftw_complex));
	zmatrix0_type.Commit();

	 ztensor0_type = zmatrix0_type.Create_hvector(n/size, 1, 										n*sizeof(fftw_complex));
	ztensor0_type.Commit();

}

void Mymat::typefree()
{
	status = 0;
	byte_type.Free();
	tensor1_type.Free();
	xtensor0_type.Free();
	ycolumn0_type.Free();
	ymatrix0_type.Free();
	ytensor0_type.Free();
	zcolumn0_type.Free();
	zmatrix0_type.Free();
	ztensor0_type.Free();
}


/* --------------------------------------------------------------------------*/
/**
* @brief 从x方向第i个进程来的数据的接受点
*
* @param i
*
* @returns   
*/
/* ----------------------------------------------------------------------------*/
double* Mymat::inbuff_x(int i)
{
	int j =i;
	j = (j+myorder[2])%size;
	return &(ele[in_position[j]][0]);
}

double* Mymat::inbuff_y(int i)
{
	int j =i;
	j = (j+myorder[1])%size;
	return &(ele[in_position[j]][0]);
}


double* Mymat::inbuff_z(int i)
{
	int j =i;
	j = (j+myorder[0])%size;
	return &(ele[in_position[j]][0]);
}


/* --------------------------------------------------------------------------*/
/**
* @brief 得到x方向上第i个进程(即xorder[i])将得到的数据的起始点
*
* @param i
*
* @returns   
*/
/* ----------------------------------------------------------------------------*/
double* Mymat::outbuff_x(int i)
{
	int j =i;
	j = (j+myorder[2])%size;
	return &(ele[out_position_x[j]][0]);
}


double* Mymat::outbuff_y(int i)
{
	int j =i;
	j = (j+myorder[1])%size;
	return &(ele[out_position_y[j]][0]);
}


double* Mymat::outbuff_z(int i)
{
	int j =i;
	j = (j+myorder[0])%size;
	return &(ele[out_position_z[j]][0]);
}


/* --------------------------------------------------------------------------*/
/**
* @brief 将各块数据传输到对应进程
*
* @param mat1
*/
/* ----------------------------------------------------------------------------*/
void Mymat::trans_x(Mymat &mat1)
{
	for(int i=0;i<size;i++)
    {
//		if(myid<=7)
//		{	
//			std::cout << myid <<"send to" << xorder[i] << " " << 
//			"receive from" << xorder[(size-i)%size] << std::endl;
//		}
    	MPI::COMM_WORLD.Sendrecv(
    			outbuff_x(i), 1, xtensor0_type, xorder[i], 99,								mat1.inbuff_x((size-i)%size), 1, tensor1_type, xorder[(size-i)%size],  99);
//		if(myid<=7)
//		{	
//			std::cout << myid <<"finish send to" << xorder[i] << " " << 
//			"receive from" << xorder[(size-i)%size] << std::endl;
//		}
    }
}

void Mymat::trans_y(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
		MPI::COMM_WORLD.Sendrecv(
			outbuff_y(i), 1, ytensor0_type, yorder[i], 99,								mat1.inbuff_y((size-i)%size), 1, tensor1_type, yorder[(size-i)%size], 99);
	}
}

void Mymat::trans_z(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
		MPI::COMM_WORLD.Sendrecv(
			outbuff_z(i), 1, ztensor0_type, zorder[i], 99,								mat1.inbuff_z((size-i)%size), 1, tensor1_type, zorder[(size-i)%size], 99);
	}
}

/* --------------------------------------------------------------------------*/
/**
* @brief 将mat1中数据按x方向返回各进程对应位置
*
* @param mat1
*/
/* ----------------------------------------------------------------------------*/
void Mymat::retrans_x(Mymat &mat1)
{
	for(int i=0;i<size;i++)
    {
    	MPI::COMM_WORLD.Sendrecv(														mat1.inbuff_x(i), 1, tensor1_type, xorder[i],  99,							outbuff_x((size-i)%size), 1, xtensor0_type, xorder[(size-i)%size], 99);
    }
}

void Mymat::retrans_y(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
    	MPI::COMM_WORLD.Sendrecv(														mat1.inbuff_y(i), 1, tensor1_type, yorder[i], 99,							outbuff_y((size-i)%size), 1, ytensor0_type, yorder[(size-i)%size], 99);
	}
}

void Mymat::retrans_z(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
    	MPI::COMM_WORLD.Sendrecv(														mat1.inbuff_z(i), 1, tensor1_type, zorder[i], 99,							outbuff_z((size-i)%size), 1, ztensor0_type, zorder[(size-i)%size], 99);
	}
}

/* --------------------------------------------------------------------------*/
/**
* @brief 沿x方向对每个向量做fft
*/
/* ----------------------------------------------------------------------------*/
void Mymat::fft(void)
{
	fftw_plan p;
	for(int i=0;i<size_m*size_n;i++)
	{
		p = fftw_plan_dft_1d(size_l, &ele[i*size_l], 											&ele[i*size_l], FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
	}
	fftw_free(p);
}


/* --------------------------------------------------------------------------*/
/**
* @brief x方向每个向量ifft
*/
/* ----------------------------------------------------------------------------*/
void Mymat::ifft(void)
{
	fftw_plan p;
	for(int i=0;i<size_m*size_n;i++)
	{
		p = fftw_plan_dft_1d(size_l, &ele[i*size_l], 											&ele[i*size_l], FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(p);
	}
	fftw_free(p);
}



/* --------------------------------------------------------------------------*/
/**
* @brief 求最大模的平方
*
* @returns l_inf^2   
*/
/* ----------------------------------------------------------------------------*/
double Mymat::norm_inf(void)
{
	double num = 0;
	double num1;
	double* p = &ele[0][0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
	{
		num1 = *p*(*p);
		p++;
		num1 += *p*(*p);
		p++;
		num = std::max(num, num1);
//		num = std::max(num, ele[i][0]*ele[i][0]+ele[i][1]*ele[i][1]);
	}
	num1 = *p*(*p);
	p++;
	num1 += *p*(*p);
	num = std::max(num, num1);
	return num;
}



Mymat& Mymat::operator+=(const Mymat& mat1)
{
	double* p  = &ele[0][0];
	double* p1 = &mat1.ele[0][0];
	for(int i=0;i<2*size_l*size_m*size_n-1;i++)
	{
		*p += *p1;
		p++;
		p1++;
//		ele[i][0] += mat1.ele[i][0];
//		ele[i][1] += mat1.ele[i][0];
	}
	*p += *p1;
	return *this;
}

Mymat& Mymat::operator/=(double alpha)
{
	double* p = &ele[0][0];
	for(int i=0;i<2*size_l*size_m*size_n-1;i++)
	{
		*p /= alpha;
		p++;
//		ele[i][0] /= alpha;
//		ele[i][1] /= alpha;
	}
	*p /= alpha;
	return *this;
}

Mymat& Mymat::operator=(const Mymat& mat1)
{
	double* p  = &ele[0][0];
	double* p1 = &mat1.ele[0][0];
	for(int i=0;i<2*size_l*size_m*size_n-1;i++)
	{
		*p = *p1;
		p++;
		p1++;
//		ele[i][0] = mat1.ele[i][0];
//		ele[i][1] = mat1.ele[i][1];
	}
	*p = *p1;
	return *this;
}

/* --------------------------------------------------------------------------*/
/**
* @brief 每个元素都自乘上其模的p次方
*
* @param p 
*
* @returns   
*/
/* ----------------------------------------------------------------------------*/
Mymat& Mymat::operator^=(int q)
{
	double val;
	double* p = &ele[0][0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
	{
		val = *(p+1);
		val = val*val;
		val += *p*(*p);
		*p *= val;
		p++;
		*p *= val;
		p++;
//		val = ele[i][0]*ele[i][0] + ele[i][1]*ele[i][1];
//		ele[i][0] *= val;
//		ele[i][1] *= val;
	}
	val = *(p+1);
	val = val*val;
	val += *p*(*p);
	*p *= val;
	p++;
	*p *= val;
	return *this;
}

Mymat Mymat::operator-(const Mymat& mat1) const
{
	Mymat mat2(size_l,size_m,size_n);
	double* p  = &ele[0][0];
	double* p1 = &mat1.ele[0][0];
	double* p2 = &mat2.ele[0][0];
	for(int i=0;i<2*size_l*size_m*size_n-1;i++)
	{
		*p2 = *p - *p1;
		p++;
		p1++;
		p2++;
//		mat2.ele[i][0] = ele[i][0] - mat1.ele[i][0];
//		mat2.ele[i][1] = ele[i][1] - mat1.ele[i][1];
	}
	*p2 = *p - *p1;
	return mat2;
}


/* --------------------------------------------------------------------------*/
/**
* @brief 为什么乘法操作通过指针实现和指标实现的结果有差异?
*
* @param alpha
*
* @returns   
*/
/* ----------------------------------------------------------------------------*/
Mymat Mymat::operator*(double alpha) const
{	
	Mymat mat2(size_l,size_m,size_n);
	double* p  = &ele[0][0];
	double* p2 = &mat2.ele[0][0];
	for(int i=0;i<2*size_l*size_m*size_n-1;i++)
	{
		*p2 = alpha * (*p);
		p2++;
		p++;
//		mat2.ele[i][0] = ele[i][0] * alpha;  
//		mat2.ele[i][1] = ele[i][1] * alpha;
	}
	*p2 = alpha * *p;
	return mat2;
}



