/*
 *  resources.h
 *  Percolation
 *
 *  Created by Samuel Schoenholz on 6/2/10.
 *  Copyright 2010 Swarthmore College. All rights reserved.
 *
 */

#ifndef RESOURCES
#define RESOURCES


#include "std_include.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <list>

//#include "Eigen/Core"

#define PI M_PI
//#define PI 3.1415926535897932384626433832795028841971693993751058


using namespace std;


template <typename T> inline T POW2(T x) {	return x*x;}
template <typename T> inline T POW3(T x) {	return x*x*x;}
template <typename T> inline T POW4(T x) 
{
	T temp = x*x;
	return temp*temp;
}
template <typename T> inline T POW6(T x) 
{
	T temp = x*x*x;
	return temp*temp;
}
template <typename T> inline T POW12(T x) 
{
	T temp = x*x*x;
	temp = temp*temp;
	return temp*temp;
}


template <typename T>
inline T max_abs_element(int N, T const *const a)
{
	T m = fabs(a[0]);
	for(int ii=1; ii<N; ii++)   if(fabs(a[ii])>m)   m=fabs(a[ii]);
	return m;
}













inline bool FileExists(std::string strFilename)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
	
	// Attempt to get the file attributes 
	intStat = stat(strFilename.c_str(),&stFileInfo);
	if(intStat == 0) {
		blnReturn = true;
	} else {
		blnReturn = false;
	}
	return(blnReturn);
}










//function to output a vector<T> as a mathematica vector
template<class T>
ostream &operator<<(ostream &out, vector<T> &data)
{
	out << "{";
	typename vector<T>::iterator it;
	for(it = data.begin();it!=data.end();it++)
		out << (*it) << ((it+1!=data.end())?"," : "");
	
	out << "}";
	return out;
}

//function to input a vector<T> as a mathematica vector
template<class T>
istream &operator>>(istream &in, vector<T> &data)
{
    data.clear();
	char current;
	in >> current;
	if(data.size()==0){
		while(current!='}')
		{
			T temp;
			in >> temp;
			data.push_back(temp);
			in >> current;
		}
	}else {
		int t = 0;
		while(current!='}')
		{
			in >> data[t]; 
			t++;
			in >> current;
		}
	}
    
	return in;
}


//function to perform arithmatic on vector<T>.
template<class T>
vector<T> &operator/(vector<T> &vec, T &div)
{
	typename vector<T>::iterator it;
	for(it = vec.begin();it!=vec.end();it++)
		(*it) = (*it)/div;
	return vec;
}

template<class T>
vector<T> &operator+(vector<T> &vec, vector<T> &add)
{
	vector<T> ret;
	for(int it = 0;it<vec.size();it++)
		ret.push_back(vec[it]+add[it]);
	return ret;
}

template<class T>
vector<T> &operator/=(vector<T> &vec,const T &div)
{
	typename vector<T>::iterator it;
	for(it = vec.begin();it!=vec.end();it++)
		(*it) = (*it)/div;
	return vec;
}

//code to do linear regressions
struct LinearRegression {
    double slope;
    double intercept;
    
    LinearRegression() : slope(0), intercept(0) {}
    LinearRegression(double _i, double _s) : slope(_s), intercept(_i) {}
};

LinearRegression LinearFit(vector<double> x_data, vector<double> y_data);

//code to compute factorials as well as absolute values and sgns.
long int factorial(long int n);

double logfactorial(int n);

double logchoose(int n, int k);

double choose( int n,  int k);

int pow(int n, int e);

template <class T>
T sgn(T a)
{
    return (a>0) - (a<0);
}

template <class T>
T abs(T a)
{
    return (a>0) ? a : -a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////CODE TO APPROXIMATE GAUSSIAN INVERSE CDF////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////From: http://www.johndcook.com/normal_cdf_inverse.html /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// compute log(1+x) without losing precision for small values of x
double LogOnePlusX(double x);

double RationalApproximation(double t);

double NormalCDFInverse(double p);

vector<string> SplitString(const string &target, const string &token);

#endif
