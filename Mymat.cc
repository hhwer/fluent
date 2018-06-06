/*
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
	plan_on = 0;
}

Mymat::Mymat(int l, int m, int n )
{
	size_l = l;
	size_m = m;
	size_n = n;
	status = 0;
	plan_on = 0;
	ele = (double*) fftw_malloc( sizeof(double)*l*m*n );
	double* p=&ele[0];
	for(int i=0;i<l*m*n-1;i++)
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
	plan_on = 0;
	double num0 = (double)num; 
	ele = (double*) malloc( sizeof(double)*l*m*n );
	double* p=&ele[0];
	for(int i=0;i<l*m*n-1;i++)
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
	plan_on = 0;
	ele = (double*) malloc( sizeof(double)*size_l*size_m*size_n );
	double* p=&ele[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
	{
		*p = 0.0;
		p++;
	}
	*p = 0.0;
}
Mymat::~Mymat(void)
	{
		fftw_free(ele);
		std::vector<double>().swap(factor);
		std::vector<int>().swap(xorder);
		std::vector<int>().swap(yorder);
		std::vector<int>().swap(zorder);
		std::vector<int>().swap(in_position);
		std::vector<int>().swap(out_position_x);
		std::vector<int>().swap(out_position_y);
		std::vector<int>().swap(out_position_z);
		if(status)
		{
			typefree();
		}
		if(plan_on)
		{
			fftw_destroy_plan(p_dct);
			fftw_destroy_plan(p_idct);
			fftw_destroy_plan(p_dst);
			fftw_destroy_plan(p_idst);
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
	byte_type = MPI::DOUBLE.Create_vector(1, 1, 1);
	byte_type.Commit();

	tensor1_type = byte_type.Create_vector(n*n/size, n, n*size);
	tensor1_type.Commit();

	xtensor0_type = byte_type.Create_vector(1, n*n*n/size, 0);
	xtensor0_type.Commit();


	ycolumn0_type = byte_type.Create_vector(n, 1, n);
	ycolumn0_type.Commit();

	ymatrix0_type = ycolumn0_type.Create_hvector(n, 1, sizeof(double));
	ymatrix0_type.Commit();

	ytensor0_type = ymatrix0_type.Create_vector(n/size, 1, 1);
	ytensor0_type.Commit();
	
	
	zcolumn0_type = byte_type.Create_vector(n, 1, n*n);
	zcolumn0_type.Commit();

	zmatrix0_type = zcolumn0_type.Create_hvector(n, 1, 											sizeof(double));
	zmatrix0_type.Commit();

	ztensor0_type = zmatrix0_type.Create_hvector(n/size, 1, 										n*sizeof(double));
	ztensor0_type.Commit();

	utensor2_type = byte_type.Create_vector(1, n*n*n, 0);
	utensor2_type.Commit();

	vcolumn2_type = byte_type.Create_vector(n, 1, n*n);
	vcolumn2_type.Commit();

	vmatrix2_type = vcolumn2_type.Create_hvector(n, 1, 											sizeof(double));
	vmatrix2_type.Commit();
	  
	vtensor2_type = vmatrix2_type.Create_hvector(n, 1, 											n*sizeof(double));
	vtensor2_type.Commit();
	
	wcolumn2_type = byte_type.Create_vector(n, 1, n);
	wcolumn2_type.Commit();

	wmatrix2_type = wcolumn2_type.Create_hvector(n, 1, 											n*n*sizeof(double));
	wmatrix2_type.Commit();
	
	wtensor2_type = wmatrix2_type.Create_hvector(n, 1, 											sizeof(double));
	wtensor2_type.Commit();
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
	
	utensor2_type.Free();
	vcolumn2_type.Free();
	vmatrix2_type.Free();
	vtensor2_type.Free();
	wcolumn2_type.Free();
	wmatrix2_type.Free();
	wtensor2_type.Free();
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
	return &(ele[in_position[j]]);
}

double* Mymat::inbuff_y(int i)
{
	int j =i;
	j = (j+myorder[1])%size;
	return &(ele[in_position[j]]);
}


double* Mymat::inbuff_z(int i)
{
	int j =i;
	j = (j+myorder[0])%size;
	return &(ele[in_position[j]]);
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
	return &(ele[out_position_x[j]]);
}


double* Mymat::outbuff_y(int i)
{
	int j =i;
	j = (j+myorder[1])%size;
	return &(ele[out_position_y[j]]);
}


double* Mymat::outbuff_z(int i)
{
	int j =i;
	j = (j+myorder[0])%size;
	return &(ele[out_position_z[j]]);
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
* @brief 从相应的块U里获得V,W
*
* @param mat1 存储V
* @param mat2 存储W
*/
/* ----------------------------------------------------------------------------*/
void Mymat::getVW(Mymat &mat1, Mymat &mat2)
{
//	std::cout<<myid<<' '<<myorder[0]*size*size+myorder[1]*size+myorder[2]<<' '<<myorder[1]*size*size+myorder[2]*size+myorder[0]<<std::endl;
    MPI::COMM_WORLD.Sendrecv(														ele, 1, vtensor2_type, 														myorder[1]*size*size+myorder[2]*size+myorder[0], 99,						mat1.ele, 1, utensor2_type, 												myorder[2]*size*size+myorder[0]*size+myorder[1], 99);
    
	MPI::COMM_WORLD.Sendrecv(														ele, 1, wtensor2_type, 														myorder[2]*size*size+myorder[0]*size+myorder[1], 99,						mat2.ele, 1, utensor2_type, 												myorder[1]*size*size+myorder[2]*size+myorder[0], 99);
}


void Mymat::getplan()
{
	plan_on = 1;
	p_dct = fftw_plan_r2r_1d(size_l, ele, 											ele, FFTW_REDFT10, FFTW_MEASURE);
	p_idct = fftw_plan_r2r_1d(size_l, ele, 											ele, FFTW_REDFT01, FFTW_MEASURE);
	p_dst = fftw_plan_r2r_1d(size_l, ele, 											ele, FFTW_RODFT10, FFTW_MEASURE);
	p_idst = fftw_plan_r2r_1d(size_l, ele, 											ele, FFTW_RODFT01, FFTW_MEASURE);
}

/* --------------------------------------------------------------------------*/
/**
* @brief 沿x方向对每个向量做fft
* 		 k=-1,1对应反对称,对称两种情形
*/
/* ----------------------------------------------------------------------------*/
void Mymat::fft(int k)
{
	if(k==1)
	{
		for(int i=0;i<size_m*size_n;i++)
		{
			fftw_execute_r2r(p_dct,&ele[i*size_l],&ele[i*size_l]);
		}
	}
	else
	{
		for(int i=0;i<size_m*size_n;i++)
		{
			fftw_execute_r2r(p_dst,&ele[i*size_l],&ele[i*size_l]);
		}
	}
	
}


/* --------------------------------------------------------------------------*/
/**
* @brief x方向每个向量ifft
* 		 k=-1,1对应反对称,对称两种情形
*/
/* ----------------------------------------------------------------------------*/
void Mymat::ifft(int k)
{
	if(k==1)
	{	
		for(int i=0;i<size_m*size_n;i++)
		{
			fftw_execute_r2r(p_idct,&ele[i*size_l],&ele[i*size_l]);
		}
	}
	else
	{
		for(int i=0;i<size_m*size_n;i++)
		{
			fftw_execute_r2r(p_idst,&ele[i*size_l],&ele[i*size_l]);
		}
	}
	
}



/* --------------------------------------------------------------------------*/
/**
* @brief 求最大模的平方的位置
*
* @returns l_inf   
*/
/* ----------------------------------------------------------------------------*/
double Mymat::norm_inf(void)
{
	std::vector <double> v;
	v.resize(size_l*size_m*size_n);
	for(int i=0;i<size_l*size_m*size_n-1;i++)
	{
		v[i] = fabs(ele[i]);
	}
	int num = (int) (std::max_element(v.begin(),v.end()) - v.begin());
	std::vector<double>().swap(v);
	return num;
}



Mymat& Mymat::operator+=(const Mymat& mat1)
{
	double* p  = &ele[0];
	double* p1 = &mat1.ele[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
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
	double* p = &ele[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
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
	double* p  = &ele[0];
	double* p1 = &mat1.ele[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
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


Mymat Mymat::operator-(const Mymat& mat1) const
{
	Mymat mat2(size_l,size_m,size_n);
	double* p  = &ele[0];
	double* p1 = &mat1.ele[0];
	double* p2 = &mat2.ele[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
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

Mymat Mymat::operator+(const Mymat& mat1) const
{
	Mymat mat2(size_l,size_m,size_n);
	double* p  = &ele[0];
	double* p1 = &mat1.ele[0];
	double* p2 = &mat2.ele[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
	{
		*p2 = *p + *p1;
		p++;
		p1++;
		p2++;
	}
	*p2 = *p + *p1;
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
	double* p  = &ele[0];
	double* p2 = &mat2.ele[0];
	for(int i=0;i<size_l*size_m*size_n-1;i++)
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


/* --------------------------------------------------------------------------*/
/**
* @brief 将某个某个数据块写入文件./b.txt
*
* @param a  进程
* @param b	文件名 
*/
/* ----------------------------------------------------------------------------*/
void Mymat::myprint(int a,int b)
{
	if(myid==a)
	{	
		std::ofstream fp;
		std::string name = "./" + std::to_string(b) + ".txt";
		fp.open(name.c_str());
		for(int k=0;k<size_n;k++)
		{	
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					fp << ele[i+j*size_l+k*size_l*size_m] << ' ';  
				}
				fp << std::endl;
				
			}
		}
		fp.close();
	}
}


