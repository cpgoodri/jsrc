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

/** @file */ 


// @file A file for commonly used global functions.


template <typename T> inline T POW2(T x){	return x*x;}	//!<x^2
template <typename T> inline T POW3(T x){	return x*x*x;}	//!<x^3
template <typename T> inline T POW4(T x) //!<x^4 
{
	T temp = x*x;
	return temp*temp;
}
template <typename T> inline T POW6(T x) //!<x^6
{
	T temp = x*x*x;
	return temp*temp;
}
template <typename T> inline T POW12(T x) //!<x^12
{
	T temp = x*x*x;
	temp = temp*temp;
	return temp*temp;
}


//! Get the largest absolute value from an array
template <typename T>
inline T max_abs_element(int N, T const *const a)
{
	T m = fabs(a[0]);
	for(int ii=1; ii<N; ii++)   if(fabs(a[ii])>m)   m=fabs(a[ii]);
	return m;
}
template <typename T>
inline T max_abs_element(Eigen::Matrix<T,Eigen::Dynamic,1> const &a)
{
	return max_abs_element(a.size(), a.data());
}

template <class ForwardIterator>
inline dbl mean(ForwardIterator first, ForwardIterator last)
{
	dbl m = (dbl)(*first);
	int Ndata = 1;
	while(++first!=last)
	{
		++Ndata;
		m += ( ((dbl)(*first))-m )/((dbl)Ndata);
	}
	return m;
}

template <class ForwardIterator>
inline int mean_var(ForwardIterator first, ForwardIterator last, dbl &m, dbl &v)
{
	m = (dbl)(*first);
	dbl m2 = (dbl)POW2( (*first) );
	int Ndata = 1;
	while(++first!=last)
	{
		++Ndata;
		m  += ( ((dbl)(*first))-m )/((dbl)Ndata);
		m2 += ( ((dbl)POW2((*first)))-m2 )/((dbl)Ndata);
	}
	v = m2 - POW2(m);
	if( v < -1e-14)
	{
		printf("WARNING: v = %e\n", v);
		assert( false );
	}
	return Ndata;
}

template <class ForwardIterator>
inline int mean_stdev(ForwardIterator first, ForwardIterator last, dbl &m, dbl &stdev)
{
	dbl v;
	int Ndata = mean_var(first, last, m, v);
	stdev = sqrt(v);
	return Ndata;
}

template <class ForwardIterator>
inline int mean_error(ForwardIterator first, ForwardIterator last, dbl &m, dbl &e)
{
	dbl v;
	int Ndata = mean_var(first, last, m, v);
	e = sqrt(v)/sqrt( (dbl)Ndata );
	return Ndata;
}



//returns true if there exists an integer x such that x^pow == n
//the return value of x is 1 if the return value is false or if pow == 0;
inline bool IsPerfectRoot(int n, int power, int &x) 
{
	x=1;
	if(power<0) return false;
	if(power==0) return n==1;
	if(power%2==0 && n<0) return false;
	
	x = round( pow(((dbl)n),1./((dbl)power)) );
	int xpow(1);
	for(int i=0; i<power; ++i)
		xpow *= x;
	if(xpow==n)
		return true;
	x=1;
	return false;
}
inline bool IsPerfectRoot(int n, int power)
{
	int x;
	return IsPerfectRoot(n, power, x); 
}


inline std::string ConvertDblToHexString(dbl x)
{
	char cstring[256];
	sprintf(cstring, "%a", x);
	return std::string(cstring);
}

inline dbl ConvertHexStringToDbl(std::string s)
{
	dbl x;
	sscanf(s.c_str(), "%lf", &x);
	return x;
}








//! Check if a file exists.
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

static void AssertThatFileExists(std::string strFilename)
{
	if(!FileExists(strFilename))
	{
		printf("WARNING: File %s not found.\n", strFilename.c_str());
		assert(false);
	}
}





//Calculate the participation ratio of a vector of length Nvar
template <int Dim>
dbl CalculateParticipationRatio(int Nvar, dbl *vec)
{
	int n = Nvar/Dim;
	assert(Nvar == n*Dim);

	dbl num, denom, pmag2;
	num = denom = 0.; 

	for(int im=0; im<n; ++im)
	{   
		pmag2 = 0.; 
		for(int dd=0; dd<Dim; ++dd)
			pmag2 += POW2(vec[Dim*im+dd]);

		num += pmag2;
		denom += pmag2*pmag2;
	}   
	return num*num/(((dbl)n)*denom);
}



//if you have a list of size N, and you want to divide it up into m partitions, set
//	s = N / m    (where this will drop the remainder)
//Each partition will have size s or s+1
template <class ForwardIterator>
void SetPartitionBaseSize(ForwardIterator first, ForwardIterator last, uint m, uint &s)
{
	uint TotalSize = last - first;
	s = TotalSize / m;
}

template <class ForwardIterator>
void GetPartitionPointer(ForwardIterator const curr, ForwardIterator const last, ForwardIterator &next, uint s)
{
	assert(last >= curr);
	uint diff = last - curr;
	uint size;
	if(diff % s == 0)	size = s;
	else				size = s+1;
	next = curr + size;
//	if(next>last)
		printf("WARNING: diff=%i, size=%i, s=%i, diff%%s=%i\n", diff, size, s, diff%s);
	assert(next <= last);
}


class EvenIterator
{
public:
	uint TargetSize;
	uint TotalSize;
	vector<uint> sizes;
	EvenIterator(uint targetsize, uint totalsize)
		: TargetSize(targetsize), TotalSize(totalsize)
	{
		initialize_sizes();
	};

	void initialize_sizes()
	{
		uint numBins = TotalSize / TargetSize;
		uint remainder = TotalSize % TargetSize;
		sizes.assign(numBins, TargetSize);

		uint bin;
		for(uint i=0; i<remainder; ++i)
		{
			bin=i%numBins;
			assert(bin >= 0 && bin < numBins);
			++sizes[bin];
		}

		uint sum=0;
		for(vector<uint>::const_iterator it=sizes.begin(); it!=sizes.end(); ++it)
			sum += (*it);
		assert(sum == TotalSize);
	}

	uint GetNumBins()
	{
		return sizes.size();
	}

	void Print()
	{
		printf("%i bins with sizes:\n", (int)sizes.size());
		for(vector<uint>::const_iterator it=sizes.begin(); it!=sizes.end(); ++it)
			printf("%i\n", (*it));

	}

};





/*
int main()
{
	vector<int> a;
	for(int i=0; i<96; ++i)	a.push_back(i);
	uint s, m = 10;
	vector<int>::iterator curr, next, last;
	curr = a.begin();
	last = a.end();
	SetPartitionBaseSize(curr, last, m, s);
	while(curr != last)
	{
		GetPartitionPointer(curr, last, next, s);
		for(vector<int>::iterator it = curr; it!=next; ++it)
			printf("%4i ", (*it));
		printf("\n");
		curr = next;
	}
}
*/






















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
    dbl slope;
    dbl intercept;
    
    LinearRegression() : slope(0), intercept(0) {}
    LinearRegression(dbl _i, dbl _s) : slope(_s), intercept(_i) {}
};

LinearRegression LinearFit(vector<dbl> x_data, vector<dbl> y_data);

//code to compute factorials as well as absolute values and sgns.
long int factorial(long int n);

dbl logfactorial(int n);

dbl logchoose(int n, int k);

dbl choose( int n,  int k);

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


//code to calculate the volume and surface area of an n-dimensional unit sphere
static dbl nSphere_Sn(int n);
static dbl nSphere_Vn(int n)
{
	if(n==0) return 1.;
	return nSphere_Sn(n-1)/((dbl)n);
}
dbl nSphere_Sn(int n)
{
	if(n==0) return 2.;
	return 2.*M_PI*nSphere_Vn(n-1);
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////CODE TO APPROXIMATE GAUSSIAN INVERSE CDF////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////From: http://www.johndcook.com/normal_cdf_inverse.html /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// compute log(1+x) without losing precision for small values of x
dbl LogOnePlusX(dbl x);

dbl RationalApproximation(dbl t);

dbl NormalCDFInverse(dbl p);

vector<string> SplitString(const string &target, char ctoken);
vector<string> SplitString(const string &target, const string &token);

#endif
