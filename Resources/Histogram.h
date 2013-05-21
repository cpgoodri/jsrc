#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "Resources.h"


template <typename T=dbl, typename U=T>
class Histogram
{
private:
	vector< pair<U,U> > hist;
	int Ndata;

	vector<T> temp;
	vector<T> const& GetD(vector<T> const &data);
	void FixedWidth(vector<T> const &data, int Nbins, U width, U min, U max);
public:
	void FixedWidth(vector<T> const &data, U width);
	void FixedWidth(vector<T> const &data, U width, U min, U max);
	void FixedWidth(vector<T> const &data, int Nbins);
	void FixedWidth(vector<T> const &data, int Nbins, U min, U max);

	void VariableWidth(vector<T> const &data, int Nppb);
	
	void MinVariableWidth(vector<T> const &data, int Nppb, U minWidth);

	//Retreive histogram:
	int Nbins() const;
	U x(int i) const;
	U y(int i) const;
	int c(int i) const;
};

template <typename T, typename U>
vector<T> const& Histogram<T,U>::GetD(vector<T> const &data)
{
	bool is_sorted = true;
	for(typename vector<T>::const_iterator it = data.begin(); it+1 != data.end(); ++it)
		if((*it) > (*(it+1)))
		{
			is_sorted = false;
			break;
		}
	if(is_sorted)
		return &data;

	temp = data;
	std::sort(temp.begin(), temp.end());
	return &temp;
}

template <typename T, typename U>
void Histogram<T,U>::FixedWdith(vector<T> const &data, int Nbins, U width, U min, U max)
{
	assert( round(((U)(max-min))/width) == Nbins );
	hist.clear();
	Ndata = (int)D.size();

	vector<T> const &D = data;//GetD(data);
	U half_width = width/2.;
	long count[Nbins];
	for(int i=0; i<Nbins; ++i)
		count[i] = 0;

	for(int i=0; i<Nbins; ++i)
		hist.push_back( std::make_pair( ((U)min) + ((U)i)*width + half_width, 0.) );
	assert((int)hist.size() == Nbins);

	int bin;
	//for(int i=0; i<(int)D.size(); ++i)
	for(typename vector<T>::const_iterator it = D.begin(); it!=D.end(); ++it)
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
void Histogram<T,U>::FixedWdith(vector<T> const &data, U width)
{
	U min = *std::min_element(data.begin(), data.end());
	U max = *std::max_element(data.begin(), data.end());
	FixedWidth(data,width,min,max);
}

template <typename T, typename U>
void Histogram<T,U>::FixedWdith(vector<T> const &data, U width, U min, U max)
{
	int Nbins = round((max-min)/width);
	FixedWidth(data,Nbins,width,min,max);
}

template <typename T, typename U>
void Histogram<T,U>::FixedWdith(vector<T> const &data, int Nbins)
{
	U min = *std::min_element(data.begin(), data.end());
	U max = *std::max_element(data.begin(), data.end());
	FixedWidth(data,Nbins,min,max);
}

template <typename T, typename U>
void Histogram<T,U>::FixedWdith(vector<T> const &data, int Nbins, U min, U max)
{
	dbl width = (max-min)/((U)Nbins);
	FixedWidth(data,Nbins,width,min,max);
}

template <typename T, typename U>
void Histogram<T,U>::VariableWidth(vector<T> const &data, int Nppb)
{
	vector<T> const &D = GetD(data);
	hist.clear();
	Ndata = (int)D.size();

	vector<T>::const_iterator it1, it2;
	int Npib; //Number of points in the current bin;
	U width, avg;
	for(int first = 0; first < (int)D.size(); first += Nppb)
	{
		Npib = std::min(Nppb,(int)D.size() - first);
		if(Npib <= 1) continue;

		it1 = D.begin()+first;
		it2 = it1 + Npib;

		width = (*it2) - (*it1);
		avg = mean(it1,it2);
		hist.push_back( std::make_pair(avg, Npib/(width*Ndata)) );
	}

	while(last_data<num_data)
	{
		data_in_bin.assign(data.begin()+first_data, data.begin()+last_data);
		width = data_in_bin.back() - data_in_bin.front();
		calc_avg(data_in_bin, avg);
		histogram.push_back( std::make_pair(avg, num_data_in_bin/(width*num_data)) );

		first_data += M;
		last_data = first_data+M-1;

}











#endif //HISTOGRAM_H
