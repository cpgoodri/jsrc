#ifndef BOND_LIST_H
#define BOND_LIST_H

#include "../Resources/std_include.h"
#include "../Resources/index_map.h"
#include <Eigen/Sparse>
#include "Bond.h"
#include "Neighbors.h"

namespace LiuJamming
{

//! Class to store a list of bonds.

/*!
 * This class stores a list of bonds and provides numerous functions for making static calculation,
 * such as computing the energy, gradient and the hessian matrix. It also provides functions for
 * removing rattlers.
 */
template<int Dim>
class CBondList
{
	typedef Eigen::Matrix<dbl, Dim, 1> dvec;
	typedef Eigen::Matrix<dbl, Dim, Dim> dmat;
	typedef CBond<Dim> BOND;
	typedef Eigen::Triplet<cdbl> TRIP;

//! @name Storage Variables
///@{
	int N;						//!<Number of possible nodes.
	dbl Volume;					//!<Volume containing the bonds.
	vector<BOND> list;		//!<A list of the bonds.

///@}

public:
//! @name Constructors and Operators
///@{
	CBondList(int _N=0, dbl V=0.);							//!<Construct a CBondsList.
	CBondList(const CBondList &src);						//!<Copy constructor.
	CBondList<Dim>& operator=(const CBondList<Dim> &src);	//!<Overloaded equals operator.

///@}

//! @name Get and Set Methods
///@{
	void AddBond(const BOND &b);	//!<Add a bond to the list
	void SetVolume(dbl _V);			//!<Set the volume.
	void SetN(int _N);				//!<Set the number of nodes.

	int GetN() const;				//!<Get the number of nodes.
	int GetNBonds() const;			//!<Get the number of bonds.
	dbl GetVolume() const;			//!<Get the volume.

///@}

//! @name Bond Manipulation
///@{
	void CalculateNeighbors(vector< vector<int> > &nbrs) const;	//!<Get a list of the neighbors of each particle
	void RemoveBonds(vector<bool> const &BondsToRemove); //!<Remove bonds from the list
	void UpdateBondIndices(index_map const &map);	//!<When some nodes are removed, as expressed by the index_map map, decrease the i and j indices of all bonds accordingly.
	int  IdentifyRattlers(vector< vector<int> > &nbrs, vector<bool> &rattlers, vector<bool> const &fixed, int c=Dim+1, bool Verbose=false) const; //!<Identify nodes that are not fixed and are involved in less than c bonds.
	int  IdentifyRattlers(vector< vector<int> > &nbrs, vector<bool> &rattlers, int c=Dim+1, bool Verbose=false) const; //!<Identify rattlers assuming no fixed nodes.
	void RemoveRattlers(int c=Dim+1, bool Verbose=false); //!<Remove rattlers, i.e. nodes with less than c bonds
	void RemoveRattlers(index_map &map, int c=Dim+1, bool Verbose=false); //!< Remove rattlers, and return the corresponding index_map.

///@}

//! @name Computations
///@{
	dbl  ComputeEnergy() const;							//!<Compute the energy.
	dbl  ComputePressure() const;						//!<Compute the pressure.
	void ComputeStressTensor(dmat &stress) const;		//!<Compute the Dim by Dim stress tensor.
//	void ComputeData(CSimpleData &data) const;
	dbl  ComputeGradient(Eigen::VectorXd &grad) const;	//!<Compute the Dim*N dimensional energy gradient (i.e. -Fnet), and return the energy.
	void ComputeHessianElements(vector<TRIP> &coefficients, dvec k, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;	//!<Compute the elements of the hessian as a list.
	void ComputeHessian(Eigen::SparseMatrix<dbl> &hess, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;					//!<Compute the hessian, note: the hessian is NOT mass-normalized.
	void ComputeHessian_BZ(Eigen::SparseMatrix<cdbl> &hess, dvec k, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;		//!<Compute the hessian at non-zero wavevector k, note: the hessian is NOT mass-normalized.

///@}

//! @name Misc.
///@{
	bool CheckConsistency() const;	//!<Check that the i and j index of every bond is less than N and greater than or equal to 0.
	void PrintBonds() const;		//!<Print the bond list to stdout.

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
CBondList<Dim>::CBondList<Dim>(int _N, dbl V)
	:N(_N), Volume(V) 
{};

template<int Dim>
CBondList<Dim>::CBondList<Dim>(const CBondList &src) 
{
	(*this) = src; 
}

/**
 *	WARNING: this method is not yet implemented!
 */
template<int Dim>
CBondList<Dim> &CBondList<Dim>::operator=(const CBondList<Dim> &src)
{
	if(this != src)
	{
		printf("THIS NEEDS TO BE IMPLEMENTED!!!\n");
	}
	return *this;
}

/**
 *	N is automatically updated so that N>b.i and N>b.j. Asserts that b.i>=0 and b.j >= 0.
 */
template<int Dim>
inline void CBondList<Dim>::AddBond(const BOND &b)
{
	assert(b.i>=0);
	assert(b.j>=0);
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
void CBondList<Dim>::SetN(int _N)
{
	N = _N;
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
void CBondList<Dim>::CalculateNeighbors(vector< vector<int> > &nbrs) const
{
	nbrs.assign(N, vector<int>());
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		nbrs[b->i].push_back(b->j);
		nbrs[b->j].push_back(b->i);
	}
}

/**
 *	Bonds that are not removed are not altered in any way, and N does not change.
 *
 *	@param[in] BondsToRemove A constant vector<bool> of length list.size() that determines
 *	which bonds to remove. Bond i is removed only if BondsToRemove[i]==true.
 */
template<int Dim>
void CBondList<Dim>::RemoveBonds(vector<bool> const &BondsToRemove)
{
	assert(BondsToRemove.size() == list.size());
	vector<BOND> tlist;
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
	for(typename vector<BOND>::iterator b=list.begin(); b!=list.end(); ++b)
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

/**
 *	@param[in,out] nbrs A vector of vectors of ints that gives the neighbors of each node. This is updated so that rattlers are removed.
 *	@param[in,out] rattlers A vector of bools such that node i is a rattler if rattlers[i]==true. Nodes can be manually designated as a rattler by setting 
 *	rattlers[i]=true before calling this function (normally, all elements of rattlers should be initialied to false). On return, this identifies the nodes that
 *	are rattlers.
 *	@param[in] fixed A vector of bools that identifies nodes as being fixed, and thus cannot become a rattler.
 *	@param[in] c The number of bonds involving a node required for that node to not be a rattler. Default is Dim+1.
 *	@param[in] Verbose Bool to determine if number of identified rattlers, etc. should be printed to stdout.
 *	@return The total number of rattlers.
 */
template<int Dim>
int CBondList<Dim>::IdentifyRattlers(vector< vector<int> > &nbrs, vector<bool> &rattlers, vector<bool> const &fixed, int c, bool Verbose) const
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
int CBondList<Dim>::IdentifyRattlers(vector< vector<int> > &nbrs, vector<bool> &rattlers, int c, bool Verbose) const
{
	vector<bool> fixed;
	fixed.assign(N,false);
	return IdentifyRattlers(nbrs,rattlers,fixed,c,Verbose);
}

template<int Dim>
void CBondList<Dim>::RemoveRattlers(int c, bool Verbose)
{
	index_map map;
	RemoveRattlers(map, c, Verbose);
}

/**
 *	This method calculates the nbrs list, identifies rattlers, sets map, removes any bonds that may involve rattlers, and updates the i and j indices of each bond (as well as N) to account for the removed rattlers.
 *	
 *	@param[out] map An index_map indicating which nodes are rattlers.
 *	@param[in] c The number of bonds involving a node required for that node to not be a rattler. Default is Dim+1.
 *	@param[in] Verbose Bool to determine if number of identified rattlers, etc. should be printed to stdout.
 */
template<int Dim>
void CBondList<Dim>::RemoveRattlers(index_map &map, int c, bool Verbose)
{
	//Calculate the nbrs list
	vector< vector<int> > nbrs;
	CalculateNeighbors(nbrs);
	
	//Identify the rattlers and set the map
	vector<bool> rattlers(N,false);
	IdentifyRattlers(nbrs, rattlers, c, Verbose);
	map.set_map(rattlers);
	
	//Remove the bonds that correspond to rattlers
	vector<bool> BondsToRemove(list.size(),false);
	for(int i=0; i<(int)list.size(); ++i)
		if(rattlers[list[i].i] || rattlers[list[i].j])
			BondsToRemove[i] = true;
	RemoveBonds(BondsToRemove);

	//Update the i and j indices of each bond (as well as N) in accordence with the map
	UpdateBondIndices(map);

	if(!CheckConsistency())
		abort();
}

template<int Dim>
dbl  CBondList<Dim>::ComputeEnergy() const
{
	dbl energy=0.;
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
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
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
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
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
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
void CBondList<Dim>::ComputeHessianElements(vector<TRIP> &coefficients, dvec k, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	dmat Fii, Fjj, Fij; //Fji = Fij.transpoze(); stressed block
	dmat Kii, Kjj, Kij; //Kji = Kij.transpose(); unstressed block
	dmat Bii, Bjj, Bij; //Bji = Bij.transpose(); B = unstress_coeff*K + stress_coeff*F

	dbl kdotr;
	cdbl eikdotr;
	int icorner, jcorner;
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
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


/**
 *	This is the primary method for calculating the hessian, which is returned as an Eigen::SparseMatrix<dbl> through the parameter hess.
 *
 *	@param[out] hess Eigen::SparseMatrix<dbl> representing the hessian.
 *	@param[in] unstress_coeff A dbl indicating the weight given to the unstressed component of the hessian. For most purposes, this should be set to 1 (default).
 *	@param[in] stress_coeff A dbl indicating the weight given to the stressed component of the hessian. For most purposes, this should be set to 1 (default). 
 *	Set this to 0 (and unstress_coeff to 1) to generate the hessian for the unstressed system.
 *	@param[in] tether A dbl indicating the strength of the tether. A value of tether is added to every diagonal element of the hessian. For most purposes, this should be set to 0 (default).
 */
template<int Dim>
void CBondList<Dim>::ComputeHessian(Eigen::SparseMatrix<dbl> &hess, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	//First, compute the complex matrix elements with k=0
	vector<TRIP> coeffs;
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

/**
 *	This method calculates the hessian at a non-zero wavevector, returning it as an Eigen::SparseMatrix<cdbl> through the parameter hess.
 *
 *	@param[out] hess Eigen::SparseMatrix<cdbl> representing the hessian.
 *	@param[in] k A dvec indicating the wavevector.
 *	@param[in] unstress_coeff A dbl indicating the weight given to the unstressed component of the hessian. For most purposes, this should be set to 1 (default).
 *	@param[in] stress_coeff A dbl indicating the weight given to the stressed component of the hessian. For most purposes, this should be set to 1 (default). 
 *	Set this to 0 (and unstress_coeff to 1) to generate the hessian for the unstressed system.
 *	@param[in] tether A dbl indicating the strength of the tether. A value of tether is added to every diagonal element of the hessian. For most purposes, this should be set to 0 (default).
 */
template<int Dim>
void CBondList<Dim>::ComputeHessian_BZ(Eigen::SparseMatrix<cdbl> &hess, dvec k, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	vector<TRIP> coeffs;
	ComputeHessianElements(coeffs, k, unstress_coeff, stress_coeff, tether);

	Eigen::SparseMatrix<cdbl> temp(Dim*N,Dim*N);
	hess = temp;
	//hess.setZero();
	hess.setFromTriplets(coeffs.begin(), coeffs.end());
	assert(hess.isCompressed());
}






template<int Dim>
void CBondList<Dim>::PrintBonds() const
{
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
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





}

#endif //BOND_LIST_H

