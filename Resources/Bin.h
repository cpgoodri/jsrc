#ifndef BIN_H
#define BIN_H

#include "std_include.h"

class CBin 
{
public:
	const int Nbin;
	const dbl min,max;
protected:
	dbl step;

public:
	CBin(int nbin, dbl _min, dbl _max) 
		: Nbin(nbin), min(_min), max(_max), step((max-min)/((dbl)Nbin))
	{
		printf("Creating bin with Nbin = %i, min = %e, max = %e, step = %e\n", Nbin, min, max, step);
	};  

	bool in_range(dbl val) const
	{
		if(val >= min && val < max)
			return true;
		return false;
	}

	int get_bin(dbl val) const
	{   
		int bin = floor((val-min)/step);
		if(bin < 0 || bin >= Nbin) printf("val = %e, min = %e, max = %e, step = %e, (val-min)/step = %e, bin = %i\n", val, min, max, step, (val-min)/step, bin);
		assert(bin >= 0); 
		assert(bin < Nbin);
		return bin;
	};  

	dbl get_bin_midpoint(int bin) const
	{   
		return min + (((dbl)bin)+0.5)*step;
	};  
};



















#endif //BIN_H
