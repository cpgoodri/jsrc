#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <cmath>
#include <assert.h>
using std::vector;
using std::pair;

template <typename T=double, typename U=T>
class Histogram
{
private:
	vector< pair<U,U> > hist;
	vector<int> count;
	int Ndata;

//	vector<T> temp;
//	vector<T> const& GetD(vector<T> const &data);
	bool CheckSort (vector<T> const &data);
	U    mean(typename vector<T>::const_iterator first, typename vector<T>::const_iterator last);

	void FixedWidth(vector<T> const &data, int Nbins, U width, U min, U max);
public:
	void FixedWidth(vector<T> const &data, U width);
	void FixedWidth(vector<T> const &data, U width, U min, U max);
	void FixedWidth(vector<T> const &data, int Nbins);
	void FixedWidth(vector<T> const &data, int Nbins, U min, U max);

private:
	void VariableWidth_Private(vector<T> const &data, int Nppb);
public:
	void VariableWidth(vector<T> const &data, int Nppb);

private:
	void MinVariableWidth_Private(vector<T> const &data, int Nppb, U minWidth);
public:
	void MinVariableWidth(vector<T> const &data, int Nppb, U minWidth);

private:
	void CDF_Private(vector<T> const &data, int Nppb);
public:
	void CDF(vector<T> const &data, int Nppb);

	//Retreive histogram:
	inline int Nbins() const			{	return (int)hist.size();	};
	inline U x(int i) const				{	return hist[i].first;		};
	inline U y(int i) const				{	return hist[i].second;		};
	inline int get_count(int i) const	{	return count[i];			};
};


template <typename T, typename U>
U Histogram<T,U>::mean(typename vector<T>::const_iterator first, typename vector<T>::const_iterator last)
{
	U m = (U)(*first);
	int ndata = 1;
	while(++first!=last)
	{   
		++ndata;
		m += ( ((U)(*first))-m )/((U)ndata);
	}   
	return m;
}

template <typename T, typename U>
bool Histogram<T,U>::CheckSort(vector<T> const &data)
{
	bool is_sorted = true;
	for(typename vector<T>::const_iterator it = data.begin(); it+1 != data.end(); ++it)
		if((*it) > (*(it+1)))
		{
			is_sorted = false;
			break;
		}
	return is_sorted;
}

template <typename T, typename U>
void Histogram<T,U>::FixedWidth(vector<T> const &data, int Nbins, U width, U min, U max)
{
	assert( round(((U)(max-min))/width) == Nbins );
	hist.clear();
	Ndata = (int)data.size();

	U half_width = width/2.;
	long count[Nbins];
	for(int i=0; i<Nbins; ++i)
		count[i] = 0;

	for(int i=0; i<Nbins; ++i)
		hist.push_back( std::make_pair( ((U)min) + ((U)i)*width + half_width, 0.) );
	assert((int)hist.size() == Nbins);

	int bin;
	//for(int i=0; i<(int)D.size(); ++i)
	for(typename vector<T>::const_iterator it = data.begin(); it!=data.end(); ++it)
	{
		bin = floor( ((U)((*it)-min))/width );

		//overflow or underflow gets attached to first and last bin
		if(bin < 0) bin = 0;
		if(bin >= Nbins) bin = Nbins-1;

		if(bin>=0 && bin < Nbins)
			++count[bin];
	}

	//Normalize the histogram
	for(int i=0; i<Nbins; ++i)
		hist[i].second = ((U)count[i])/(((U)data.size())*width);
}

template <typename T, typename U>
void Histogram<T,U>::FixedWidth(vector<T> const &data, U width)
{
	U min = *std::min_element(data.begin(), data.end());
	U max = *std::max_element(data.begin(), data.end());
	FixedWidth(data,width,min,max);
}

template <typename T, typename U>
void Histogram<T,U>::FixedWidth(vector<T> const &data, U width, U min, U max)
{
	int Nbins = round((max-min)/width);
	FixedWidth(data,Nbins,width,min,max);
}

template <typename T, typename U>
void Histogram<T,U>::FixedWidth(vector<T> const &data, int Nbins)
{
	U min = *std::min_element(data.begin(), data.end());
	U max = *std::max_element(data.begin(), data.end());
	FixedWidth(data,Nbins,min,max);
}

template <typename T, typename U>
void Histogram<T,U>::FixedWidth(vector<T> const &data, int Nbins, U min, U max)
{
	U width = (max-min)/((U)Nbins);
	FixedWidth(data,Nbins,width,min,max);
}

template <typename T, typename U>
void Histogram<T,U>::VariableWidth_Private(vector<T> const &data, int Nppb)
{
	assert(CheckSort(data));
	hist.clear();
	count.clear();
	Ndata = (int)data.size();

	EvenIterator ei(Nppb,Ndata);
//	ei.Print();

	uint numBins = ei.GetNumBins();
	typename vector<T>::const_iterator curr, next, last;
	curr = data.begin();
	last = data.end();
	int Npib; //Number of points in the current bin;
	U width, avg;
	for(uint bi=0; bi<numBins; ++bi)
	{
		Npib = ei.sizes[bi];
		next = curr + Npib;
		assert(next <= last);

		width = (*(next-1)) - (*curr);
		avg = mean(curr,next);
		hist.push_back( std::make_pair(avg, Npib/(width*Ndata)) );
		count.push_back( Npib );

		curr = next;
	}

/*
	int s;
	if(Ndata%Nppb==0)
		s = Nppb;
	else
		s = Nppb-1;

	typename vector<T>::const_iterator curr, next, last;
	curr = data.begin();
	last = data.end();
	int Npib; //Number of points in the current bin;
	U width, avg;
	while(curr != last)
	{
		GetPartitionPointer(curr, last, next, s);
		//current partition goes from curr to next.
		Npib = next - curr;
		assert(Npib == Nppb || Npib == Nppb-1);
	
		width = (*(next-1)) - (*curr);
		avg = mean(curr,next);
		hist.push_back( std::make_pair(avg, Npib/(width*Ndata)) );
		count.push_back( Npib );

		curr = next;
	}
*/

/*	
	typename vector<T>::const_iterator it1, it2;
	int Npib; //Number of points in the current bin;
	U width, avg;
	for(int first = 0; first < Ndata; first += Nppb)
	{
		Npib = std::min(Nppb,Ndata - first);
		if(Npib <= 1) continue;

		it1 = data.begin()+first;
		it2 = it1 + Npib; //This points to the first data element beyond the bin

		width = (*(it2-1)) - (*it1);
		avg = mean(it1,it2);
		hist.push_back( std::make_pair(avg, Npib/(width*Ndata)) );
		count.push_back( Npib );
	}
*/
	assert(hist.size() == count.size());
}

template <typename T, typename U>
void Histogram<T,U>::VariableWidth(vector<T> const &data, int Nppb)
{
	if(CheckSort(data))
		VariableWidth_Private(data, Nppb);
	else
	{
		vector<T> sdata = data;
		std::sort(sdata.begin(), sdata.end());
		VariableWidth_Private(sdata, Nppb);
	}
}

/*
template <typename T, typename U>
void Histogram<T,U>::MinVariableWidth_Private(vector<T> const &data, int Nppb, U minWidth)
{
	assert(CheckSort(data));
	hist.clear();
	Ndata = (int)data.size();

	typename vector<T>::const_iterator it1, it2;
	int Npib; //Number of points in the current bin;
	U width, avg;
	for(int first = 0; first < Ndata; )
	{
		Npib = std::min(Nppb,Ndata - first);
		if(Npib <= 1)
		{
			first += Npib;
			continue;
		}

		it1 = data.begin()+first;
		it2 = std::lower_bound(it1+Npib, data.end(), (*it1)+minWidth);
//		it2 = it1 + Npib;
		width = (*it2) - (*it1);
		Npib = (int)(it2-it1);

		avg = mean(it1,it2);
		hist.push_back( std::make_pair(avg, Npib/(width*Ndata)) );

		first += Npib;
	}
}

template <typename T, typename U>
void Histogram<T,U>::MinVariableWidth(vector<T> const &data, int Nppb, U minWidth)
{
	if(CheckSort(data))
		MinVariableWidth_Private(data, Nppb, minWidth);
	else
	{
		vector<T> sdata = data;
		std::sort(sdata.begin(), sdata.end());
		MinVariableWidth_Private(sdata, Nppb, minWidth);
	}
}
*/

template <typename T, typename U>
void Histogram<T,U>::CDF_Private(vector<T> const &data, int Nppb)
{
	assert(CheckSort(data));
	hist.clear();
	count.clear();
	Ndata = (int)data.size();

	typename vector<T>::const_iterator it;
	uint n = 0;
	for(it = data.begin(); it!=data.end(); ++it)
	{
		++n;
		if(n%Nppb==Nppb-1 || n==Ndata)
		{
			hist.push_back( std::make_pair((*it),((dbl)n)/((dbl)Ndata)) );
			count.push_back( n );
		}
	}
	assert(hist.size() == count.size());
}

template <typename T, typename U>
void Histogram<T,U>::CDF(vector<T> const &data, int Nppb)
{
	if(CheckSort(data))
		CDF_Private(data, Nppb);
	else
	{
		vector<T> sdata = data;
		std::sort(sdata.begin(), sdata.end());
		CDF_Private(sdata, Nppb);
	}
}







#endif //HISTOGRAM_H
