#ifndef BOND_LIST_H
#define BOND_LIST_H

#include "../Resources/std_include.h"
#include "../Resources/index_map.h"
#include <Eigen/Sparse>
#include "Bond.h"
#include "Neighbors.h"


template<int Dim>
class CBondList
{
	typedef Eigen::Matrix<dbl, Dim, 1> dvec;
	typedef Eigen::Matrix<dbl, Dim, Dim> dmat;
	typedef CBond<Dim> BOND;
	typedef Eigen::Triplet<cdbl> TRIP;

	int N;
	dbl Volume;
	std::vector<BOND> list;

public:
	CBondList(int _N=0, dbl V=0.)
		:N(_N), Volume(V) {};

	CBondList(const CBondList &src) { (*this) = src; };
	CBondList<Dim>& operator=(const CBondList<Dim> &src);

	void AddBond(const BOND &b);
	void SetVolume(dbl _V);

	int GetN() const;
	int GetNBonds() const;
	dbl GetVolume() const;

	//Get a list of the neighbors of each particle
	void CalculateNeighbors(std::vector< std::vector<int> > &nbrs) const;

	//Remove Bonds/Rattlers
	void RemoveBonds(std::vector<bool> const &BondsToRemove);
	void UpdateBondIndices(index_map const &map);
	int  IdentifyRattlers(std::vector< std::vector<int> > &nbrs, std::vector<bool> &rattlers, std::vector<bool> const &fixed, int c=Dim+1, bool Verbose=false) const;
	int  IdentifyRattlers(std::vector< std::vector<int> > &nbrs, std::vector<bool> &rattlers, int c=Dim+1, bool Verbose=false) const;
	void RemoveRattlers(int c=Dim+1, bool Verbose=false);
	void RemoveRattlers(index_map &map, int c=Dim+1, bool Verbose=false);

	//Compute data
	dbl  ComputeEnergy() const;
	dbl  ComputePressure() const;
	void ComputeStressTensor(dmat &stress) const;
//	void ComputeData(CSimpleData &data) const;

	//Compute gradient
	dbl  ComputeGradient(Eigen::VectorXd &grad) const; //Returns the energy

	//Compute the hessian and dynamical matrix
private:
	void ComputeHessianElements(std::vector<TRIP> &coefficients, dvec k, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;
public:
	void ComputeHessian(Eigen::SparseMatrix<dbl> &hess, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;  //hess is NOT mass-normalized
	void ComputeHessian_BZ(Eigen::SparseMatrix<cdbl> &hess, dvec k, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;  //hess is NOT mass-normalized

	bool CheckConsistency() const;
	void PrintBonds() const;
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
	if(b.i>=N) N=b.i+1;
	if(b.j>=N) N=b.j+1;
	list.push_back(b);
}

template<int Dim>
void CBondList<Dim>::SetVolume(dbl _V)
{
	Volume = _V;
}

template<int Dim>
int CBondList<Dim>::GetN() const
{
	return N;
}

template<int Dim>
int CBondList<Dim>::GetNBonds() const
{
	return (int)list.size();
}

template<int Dim>
dbl CBondList<Dim>::GetVolume() const
{
	return Volume;
}

template<int Dim>
void CBondList<Dim>::CalculateNeighbors(std::vector< std::vector<int> > &nbrs) const
{
	nbrs.assign(N, std::vector<int>());
	for(typename std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		nbrs[b->i].push_back(b->j);
		nbrs[b->j].push_back(b->i);
	}
}

//Remove a set of bonds from list. Do not change any of the bonds (like their indices)
//BondsToRemove[i] = true if the i'th bond should be removed, false otherwise
template<int Dim>
void CBondList<Dim>::RemoveBonds(std::vector<bool> const &BondsToRemove)
{
	assert(BondsToRemove.size() == list.size());
	std::vector<BOND> tlist;
	tlist.reserve(list.size());

	for(int i=0; i<(int)list.size(); ++i)
		if(!BondsToRemove[i])
			tlist.push_back(list[i]);

	tlist.swap(list);
}
	
template<int Dim>
void CBondList<Dim>::UpdateBondIndices(index_map const &map)
{
	assert(map.full_size == N);
	for(typename std::vector<BOND>::iterator b=list.begin(); b!=list.end(); ++b)
	{
		//check that the bond does not involve a rattler
		assert(map.inv(b->i) != -1);
		assert(map.inv(b->j) != -1);

		//change the i and j index according to the map
		b->i = map.inv(b->i);
		b->j = map.inv(b->j);
	}
	N = map.size();
}

//rattlers should already be initialized. if rattler[i]==true to begin, it is assumed to already be designated a rattler
//return the number of total rattlers
template<int Dim>
int CBondList<Dim>::IdentifyRattlers(std::vector< std::vector<int> > &nbrs, std::vector<bool> &rattlers, std::vector<bool> const &fixed, int c, bool Verbose) const
{
	assert(rattlers.size() == nbrs.size());
	assert(rattlers.size() == N);
	int num_total_rattlers = 0;
	int rattlers_found = 0;
	int new_rattlers_found;
	if(Verbose)
		printf("Removing particles with less than %i contacts...\t", c);

	//If a particle is already designated a rattler, clear its nbrs list...
	for(int i=0; i<N; ++i) 
		if(rattlers[i])
		{
			++num_total_rattlers;
			nbrs[i].clear();
		}

	// ... and remove it from other particle's nbrs list.
	for(int i=0; i<N; ++i) 
	{    
		for(int j=(int)nbrs[i].size()-1; j>=0; --j) 
		{    
			if(rattlers[nbrs[i][j]])
				nbrs[i].erase(nbrs[i].begin()+j);
		}    
	} 
    //Look for new rattlers recursively.
	do{
		new_rattlers_found = 0;
		for(int i=0; i<N; ++i)
		{
			if((int)nbrs[i].size() < c && !rattlers[i] && !fixed[i])
			{
				RemoveFromNeighborsList(nbrs,i);
				rattlers[i] = true;
				num_total_rattlers++;
				rattlers_found++;
				new_rattlers_found++;
			}
		}
	}while(new_rattlers_found);
	if(Verbose)
		printf("    --Removed % 5i new rattlers--\t", rattlers_found);

	//Perform Checks
	for(int i=0; i<N; ++i)
	{
		if(!fixed[i])
		{
			if(rattlers[i] && ((int)nbrs[i].size())>0)
			{
				printf("Error removing rattlers: rattlers have neighbors!\n");
				exit(EXIT_FAILURE);
			}   
			if(!rattlers[i] && ((int)nbrs[i].size())<c)
			{
				printf("Error removing rattlers: non-rattlers have too few neighbors!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	//Check the consistency of the nbrs list
	if(!CheckNeighborsConsistency(nbrs))
	{
		printf("Error: nbrs list is not consistent!\n");
		exit(EXIT_FAILURE);
	}

	if(Verbose)
		 printf("    --% 5i total rattlers--\n", num_total_rattlers);

	return num_total_rattlers;
}

template<int Dim>
int CBondList<Dim>::IdentifyRattlers(std::vector< std::vector<int> > &nbrs, std::vector<bool> &rattlers, int c, bool Verbose) const
{
	std::vector<bool> fixed;
	fixed.assign(N,false);
	return IdentifyRattlers(nbrs,rattlers,fixed,c,Verbose);
}

template<int Dim>
void CBondList<Dim>::RemoveRattlers(int c, bool Verbose)
{
	index_map map;
	RemoveRattlers(map, c, Verbose);
}

template<int Dim>
void CBondList<Dim>::RemoveRattlers(index_map &map, int c, bool Verbose)
{
	//Calculate the nbrs list
	std::vector< std::vector<int> > nbrs;
	CalculateNeighbors(nbrs);
	
	//Identify the rattlers and set the map
	std::vector<bool> rattlers(N,false);
	IdentifyRattlers(nbrs, rattlers, c, Verbose);
	map.set_map(rattlers);
	
	//Remove the bonds that correspond to rattlers
	std::vector<bool> BondsToRemove(list.size(),false);
	for(int i=0; i<(int)list.size(); ++i)
		if(rattlers[list[i].i] || rattlers[list[i].j])
			BondsToRemove[i] = true;
	RemoveBonds(BondsToRemove);

	//Update the i and j indices of each bond (as well as N) in accordence with the map
	UpdateBondIndices(map);
}

template<int Dim>
dbl  CBondList<Dim>::ComputeEnergy() const
{
//	printf("110\n"); fflush(stdout);
	dbl energy=0.;
//	for(typename std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
//		energy += b->E;
	for(int b=0; b<(int)list.size(); ++b)
	{
//		printf("111 %i\n", b); fflush(stdout);
		energy += list[b].E;
	}
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
	for(typename std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
		stress += b->g*((b->r)*(b->r.transpose()))/(b->rlen);
	stress /= Volume;
}


template<int Dim>
dbl  CBondList<Dim>::ComputeGradient(Eigen::VectorXd &grad) const
{
	//CheckConsistency();
	grad.setZero(N*Dim);
	dbl energy = 0.;
	dvec temp;
	int dd;
	for(typename std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
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

template<int Dim>
void CBondList<Dim>::ComputeHessianElements(std::vector<TRIP> &coefficients, dvec k, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	dmat Fii, Fjj, Fij; //Fji = Fij.transpoze(); stressed block
	dmat Kii, Kjj, Kij; //Kji = Kij.transpose(); unstressed block
	dmat Bii, Bjj, Bij; //Bji = Bij.transpose(); B = unstress_coeff*K + stress_coeff*F

	dbl kdotr;
	cdbl eikdotr;
	int icorner, jcorner;
	for(typename std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		//Calculate exp{i k\cdot r} = cos(k \cdot r) + i sin(k \cdot r)
		kdotr = k.dot(b->r);
		eikdotr = cdbl(cos(kdotr), sin(kdotr));

		//Calculate the matrix blocks
		b->CalculateMatrixBlocks(Fij,Kij); // equal to: (*b).CalculateMatrixBlocks(Fij,Kij);
		Fii = Fjj = -Fij;
		Kii = Kjj = -Kij;

		Bii = unstress_coeff*Kii + stress_coeff*Fii;
		Bjj = unstress_coeff*Kjj + stress_coeff*Fjj;
		Bij = unstress_coeff*Kij + stress_coeff*Fij;

		//Add the matrix blocks to the list of coefficients
		icorner = Dim*b->i;
		jcorner = Dim*b->j;
		for(int d1=0; d1<Dim; ++d1)
			for(int d2=0; d2<Dim; ++d2)
			{
				coefficients.push_back( TRIP(icorner+d1, icorner+d2, Bii(d1,d2)) );
				coefficients.push_back( TRIP(jcorner+d1, jcorner+d2, Bjj(d1,d2)) );
				coefficients.push_back( TRIP(icorner+d1, jcorner+d2, eikdotr*Bij(d1,d2)) );
				coefficients.push_back( TRIP(jcorner+d1, icorner+d2, std::conj(eikdotr)*(Bij.transpose())(d1,d2)) );
			}
	}
	
	//add the tether
	for(int ii=0; ii<Dim*N; ++ii)
		coefficients.push_back( TRIP(ii, ii, cdbl(tether,0)) );
}

template<int Dim>
void CBondList<Dim>::ComputeHessian(Eigen::SparseMatrix<dbl> &hess, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	//First, compute the complex matrix elements with k=0
	std::vector<TRIP> coeffs;
	dvec k = dvec::Zero();
	ComputeHessianElements(coeffs, k, unstress_coeff, stress_coeff, tether);
	
	//Create a temporary complex matrix
//	Eigen::SparseMatrix<cdbl> hess_temp;
	Eigen::SparseMatrix<cdbl> hess_temp(Dim*N,Dim*N);
	hess_temp.setFromTriplets(coeffs.begin(), coeffs.end());

	//Now, copy the real part to hess
	hess.setZero();
	hess = hess_temp.real();
	assert(hess.isCompressed());
}

template<int Dim>
void CBondList<Dim>::ComputeHessian_BZ(Eigen::SparseMatrix<cdbl> &hess, dvec k, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	std::vector<TRIP> coeffs;
	ComputeHessianElements(coeffs, k, unstress_coeff, stress_coeff, tether);

	hess.setZero();
	hess.setFromTriplets(coeffs.begin(), coeffs.end());
	assert(hess.isCompressed());
}






template<int Dim>
void CBondList<Dim>::PrintBonds() const
{
	for(typename std::vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
		b->print();
}

template<int Dim>
bool CBondList<Dim>::CheckConsistency() const
{
	//Calculate the largest and smallest particle index
	int imin, imax;
	imin = imax = list[0].i;
	for(int b=1; b<(int)list.size(); ++b)
	{
		if(imin > list[b].i) imin = list[b].i;
		if(imax < list[b].i) imax = list[b].i;
	}
	for(int b=0; b<(int)list.size(); ++b)
	{
		if(imin > list[b].j) imin = list[b].j;
		if(imax < list[b].j) imax = list[b].j;
	}

	printf("imin = %i, imax = %i, N = %i\n", imin, imax, N);

	if(imin >= 0 && imax < N)
		return true;
	return false;
}



//template<int Dim>
//void CBondList<Dim>::
//template<int Dim>
//void CBondList<Dim>::

#endif //BOND_LIST_H

