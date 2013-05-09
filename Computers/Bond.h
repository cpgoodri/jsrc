#ifndef BOND_H
#define BOND_H

#include "../Resources/std_include.h"


//! Class to store a single bond.

template<int Dim>
class CBond
{
	typedef Eigen::Matrix<dbl, Dim, 1> dvec;
	typedef Eigen::Matrix<dbl, Dim, Dim> dmat;
public:
	int	i;		//!<Index i.
	int j;		//!<Index j.
	dvec r;		//!<A dvec pointing from i to j.
	dbl sigma;	//!<Relevant length scale.
	dbl rlen;	//!<Distance between i and j.
	dbl E;		//!<Potential energy, V, stored in the bond.
	dbl g;		//!<First derivative of V w.r.t. rlen.
	dbl k;		//!<Second derivative of V w.r.t rlen.


//! @name Constructors and Operators
///@{
	CBond(int _i=0, int _j=0)
		:i(_i), j(_j), sigma(0.), rlen(0.), E(0.), g(0.), k(0.)
		{r = dvec::Zero(); };
	CBond(int _i, int _j, dbl _sigma, dbl _rlen, dbl _E, dbl _g, dbl _k, const dvec &_r)
		:i(_i), j(_j), sigma(_sigma), rlen(_rlen), E(_E), g(_g), k(_k)
		{r = _r; };
	CBond(const CBond &src) { *this = src; };
	CBond<Dim>& operator=(const CBond<Dim> &src);

///@}

//! @name Misc.
///@{
	void CalculateMatrixBlocks(dmat &Fij, dmat &Kij) const; //!<Calculate blocks for the hessian.
	void print() const;	//!<Print the bond to the terminal.

///@}
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

/**
 *	@param[out] Fij The stressed component of the matrix block.
 *	@param[out] Kij The unstressed component of the matrix block.
 */
template<int Dim>
void CBond<Dim>::CalculateMatrixBlocks(dmat &Fij, dmat &Kij) const
{
	Fij = -g*(dmat::Identity() - r*(r.transpose())/POW2(rlen))/rlen;
	Kij = -k*r*(r.transpose())/POW2(rlen);
}

template<int Dim>
void CBond<Dim>::print() const
{
	printf("bond between %5i and %5i: r = ", i, j);
	for(int dd=0; dd<Dim; ++dd)
		printf("% e ", r[dd]);
	printf("rlen =% e, sigma =% e, E =% e, g =% e, k =% e\n", rlen, sigma, E, g, k);
}













#endif //BOND_H

