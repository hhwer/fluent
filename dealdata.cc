
/**
* @file dealdata.cc
* @brief 数据处理
* @author hh
* @version 
* @date 2018-06-07
*/



#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>

int main(int argc, char** argv)
{

	std::ifstream infile;
	std::ofstream outfile;

	double *num,*num2;
	num = (double*)malloc(100000*sizeof(double));
	num2 = (double*)malloc(100000*sizeof(double));
	int NN=0;


	double a;
	int i =0;
	
	infile.open("normW.txt",std::ios::in);
	
	if(!infile.is_open())
		std::cout << "open infile failure" << std::endl;

	while(!infile.eof())
	{
		infile >> num[i];
		i++;
	}
	NN = i-1;

	infile.close();
//
//
	infile.open("normW1.txt",std::ios::in);
	if(!infile.is_open())
		std::cout << "open infile failure" << std::endl;

	i=0;	
	while(!infile.eof())
	{
		for(int j=0;j<9;j++){
			infile >>a;
			}
		infile >> num2[i];
		i++;
	}
	infile.close();

	outfile.open("distance.txt",std::ios::app);
	for( int i =0;i<NN;i++)
	{
		outfile << std::setprecision(32)<<num2[i]-num[i] << std::endl;		
	}
	outfile.close();





////
////
////
//	outfile.open("logW.txt",std::ios::app);
//	if(!outfile.is_open())
//		std::cout << "open file failure" << std::endl;
//	for( i =0;i<NN-1;i++)
//	{
//		outfile << std::setprecision(32)<<log(num[i+1])-log(num[i]) << std::endl;		
//	}
//	outfile.close();
//
//
//
//	outfile.open("slope3W.txt",std::ios::app);
//	if(!outfile.is_open())
//		std::cout << "open file failure" << std::endl;
//	for(int i =0;i<NN-1;i++)
//	{
//		outfile << std::setprecision(32)<< num[i+1]-num[i] << std::endl;		
//	}
//	outfile.close();
//
//
//
//	infile.open("slopeW.txt",std::ios::in);
//	
//	if(!infile.is_open())
//		std::cout << "open infile failure" << std::endl;
//
//	i =0;
//	while(!infile.eof())
//	{
//		infile >> num[i];
//		i++;
//	}
//	infile.close();
//	
//
//
//	outfile.open("slope2W.txt",std::ios::app);
//	if(!outfile.is_open())
//		std::cout << "open file failure" << std::endl;
//	for(int i =0;i<NN-2;i++)
//	{
//		outfile << std::setprecision(32)<< num[i+1]-num[i] << std::endl;		
//	}
//	outfile.close();
//
//
//
//	infile.open("logW.txt",std::ios::in);
//	
//	if(!infile.is_open())
//		std::cout << "open infile failure" << std::endl;
//
//	i =0;
//	while(!infile.eof())
//	{
//		infile >> num[i];
//		i++;
//	}
//	infile.close();
//	
//
//
//	outfile.open("log2W.txt",std::ios::app);
//	if(!outfile.is_open())
//		std::cout << "open file failure" << std::endl;
//	for(int i =0;i<NN-2;i++)
//	{
//		outfile << std::setprecision(32)<< num[i+1]-num[i] << std::endl;		
//	}
//	outfile.close();
//
//
//
//	infile.open("normU.txt",std::ios::in);	
//	if(!infile.is_open())
//		std::cout << "open infile failure" << std::endl;
//
//	i=0;
//	while(!infile.eof())
//	{
//		infile >> num[i];
//		i++;
//	}
//	infile.close();
//	
//
//
//	outfile.open("logU.txt",std::ios::app);
//	if(!outfile.is_open())
//		std::cout << "open file failure" << std::endl;
//	for(int i =0;i<NN-1;i++)
//	{
//		outfile << std::setprecision(32)<< log(num[i+1])-log(num[i]) << std::endl;		
//	}
//	outfile.close();
//	
//
//
//	outfile.open("slopeU.txt",std::ios::app);
//	if(!outfile.is_open())
//		std::cout << "open file failure" << std::endl;
//	for(int i =0;i<NN-1;i++)
//	{
//		outfile << std::setprecision(32)<< num[i+1]-num[i] << std::endl;		
//	}
//	outfile.close();
//	
//	infile.open("slopeU.txt",std::ios::in);
//	
//	if(!infile.is_open())
//		std::cout << "open infile failure" << std::endl;
//
//	i =0;
//	while(!infile.eof())
//	{
//		infile >> num[i];
//		i++;
//	}
//	infile.close();
//	
//
//
//	outfile.open("slope2U.txt",std::ios::app);
//	if(!outfile.is_open())
//		std::cout << "open file failure" << std::endl;
//	for(int i =0;i<NN-2;i++)
//	{
//		outfile << std::setprecision(32)<< num[i+1]-num[i] << std::endl;		
//	}
//	outfile.close();
//
//
//
//	infile.open("logU.txt",std::ios::in);
//	
//	if(!infile.is_open())
//		std::cout << "open infile failure" << std::endl;
//
//	i =0;
//	while(!infile.eof())
//	{
//		infile >> num[i];
//		i++;
//	}
//	infile.close();
//	
//
//
//	outfile.open("log2U.txt",std::ios::app);
//	if(!outfile.is_open())
//		std::cout << "open file failure" << std::endl;
//	for(int i =0;i<NN-2;i++)
//	{
//		outfile << std::setprecision(32)<< num[i+1]-num[i] << std::endl;		
//	}
//	outfile.close();

	free(num);
}


