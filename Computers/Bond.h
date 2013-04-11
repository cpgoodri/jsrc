#ifndef BOND_H
#define BOND_H

#include "../Resources/std_include.h"

template<int Dim>
class CBond
{
	typedef Eigen::Matrix<dbl, Dim, 1> dvec;
	typedef Eigen::Matrix<dbl, Dim, Dim> dmat;
public:
	int i,j;
	dvec r;
	dbl sigma;
	dbl rlen;
	dbl E, g, k;

	CBond(int _i=0, int _j=0)
		:i(_i), j(_j), sigma(0.), rlen(0.), E(0.), g(0.), k(0.)
		{r = dvec::Zero(); };

	CBond(int _i, int _j, dbl _sigma, dbl _rlen, dbl _E, dbl _g, dbl _k, const dvec &_r)
		:i(_i), j(_j), sigma(_sigma), rlen(_rlen), E(_E), g(_g), k(_k)
		{r = _r; };

	CBond(const CBond &src) { *this = src; };

	CBond<Dim>& operator=(const CBond<Dim> &src);

	//Print the bond to the terminal.
	void print() const;
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template<int Dim>
CBond<Dim> &CBond<Dim>::operator=(const CBond<Dim> &src)
{
	if(this != &src)
	{
		i = src.i;
		j = src.j;
		r = src.r;
		sigma = src.sigma;
		rlen = src.rlen;
		E = src.E;
		g = src.g;
		k = src.k;
	}
	return *this;
}

template<int Dim>
void CBond<Dim>::print() const
{
	printf("CBond<%i> between %5i and %5i: r = ", Dim, i, j);
	for(int dd=0; dd<Dim; ++dd)
		printf("% e ", r[dd]);
	printf("rlen =% e, sigma =% e, E =% e, g =% e, k =% e\n", rlen, sigma, E, g, k);
}













#endif //BOND_H

