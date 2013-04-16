#ifndef BOND_LIST_H
#define BOND_LIST_H

#include "../Resources/std_include.h"
#include <Eigen/Sparse>
#include "Bond.h"

template<int Dim>
class CBondList
{
	typedef Eigen::Matrix<dbl, Dim, 1> dvec;
	typedef Eigen::Matrix<dbl, Dim, Dim> dmat;
	typedef CBond<Dim> BOND;
public:

	int N;
	dbl Volume;
	std::vector<BOND> list;

	CBondList(int _N=0, dbl V=0.)
		:N(_N), Volume(V) {};

	CBondList(const CBondList &src) { (*this) = src; };
	CBondList<Dim>& operator=(const CBondList<Dim> &src);

	void AddBond(const BOND &b);

	//Compute data
	dbl  ComputeEnergy() const;
	dbl  ComputePressure() const;
	void ComputeStressTensor(dmat &stress) const;
//	void ComputeData(CSimpleData &data) const;

	//Compute gradient and hessian
	dbl  ComputeGradient(Eigen::VectorXd &grad) const; //Returns the energy
	void ComputeHessian(Eigen::SparseMatrix<dbl> &hess) const;  //hess is NOT mass-normalized
//	void ComputeHessian_BZ(Eigen::SparseMatrix< std::complex<dbl> > &hess, dvec k) const;  //hess is NOT mass-normalized
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template<int Dim>
CBondList<Dim> &CBondList<Dim>::operator=(const CBondList<Dim> &src)
{
	if(this != src)
	{
		printf("THIS NEEDS TO BE IMPLEMENTED!!!\n");
	}
	return *this;
}

template<int Dim>
void CBondList<Dim>::AddBond(const BOND &b)
{
	if(b.i>N) N=b.i;
	if(b.j>N) N=b.j;
	list.push_back(b);
}

template<int Dim>
dbl  CBondList<Dim>::ComputeEnergy() const
{
	dbl energy;
	for(std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
		energy += b->E;
	return energy;
}

template<int Dim>
dbl  CBondList<Dim>::ComputePressure() const
{
	dmat stress;
	ComputeStressTensor(stress);
	return -stress.trace()/((dbl)Dim);
}

template<int Dim>
void CBondList<Dim>::ComputeStressTensor(dmat &stress) const
{
	stress.setZero();
	for(std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
		stress += b->g*((b->r)*(b->r.transpose()))/(b->rlen);
	stress /= Volume;
}


template<int Dim>
dbl  CBondList<Dim>::ComputeGradient(Eigen::VectorXd &grad) const
{
	grad.setZero();
	dbl energy = 0.;
	dvec temp;
	int dd;
	for(std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		energy += b->E;
		temp = (b->g/b->rlen)*b->r;
		for(dd=0; dd<Dim; ++dd)
		{
			grad[Dim*b->i+dd] += temp[dd]; //Check signs here!!!
			grad[Dim*b->j+dd] -= temp[dd];
		}
	}
	return energy;
}

//template<int Dim>
//void CBondList<Dim>::
//template<int Dim>
//void CBondList<Dim>::
//template<int Dim>
//void CBondList<Dim>::
//template<int Dim>
//void CBondList<Dim>::













#endif //BOND_LIST_H

