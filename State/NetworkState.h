#ifndef NETWORK_STATE
#define NETWORK_STATE

#include "../State/StaticState.h"
#include "../Computers/StaticComputer.h"
#include "../Computers/BondList.h"
#include "../Potentials/HarmonicSpringPotential.h"



namespace LiuJamming
{

template <int Dim> class CNetworkDatabase;

template <int Dim>
class CNetworkState
{
private:
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	typedef Eigen::Matrix<int,Dim,1> ivec;

	int N;							//!<Number of nodes
	int Nbonds;						//!<Number of bonds
	Eigen::VectorXd Positions;		//!<dN vector of positions (in normalized coordinates)
	Eigen::VectorXi Bondi;
	Eigen::VectorXi Bondj;
	Eigen::VectorXd Stiffnesses;	//!<Nbonds vector of spring stiffnesses
	Eigen::VectorXd ELengths;		//!<Nbonds vector of equilibrium lengths (in real units)

	CBox<Dim> *Box;					//!<Pointer to a CBox object
	CHarmonicSpringPotential *Potential;	//!<Pointer to a CHarmonicSpringPotential object
		
	
public:
//Constructors/Destructors and copy operators
	CNetworkState();														//!<Default constructor
//	CNetworkState(int N, int Nb, dbl V);									//!<"empty" constructor
	CNetworkState(CStaticState<Dim> const &s_in, bool MakeUnstressed, bool UseUnitStiffness);
	CNetworkState(CNetworkState const &copy);								//!<Copy constructor
	CNetworkState(CNetworkState const &copy, ivec copies_per_side);			//!<Construct the CNetworkState from multiple copies of copy.
	CNetworkState<Dim> &operator=(CNetworkState<Dim> const &copy);			//!<Copy operator
	~CNetworkState();

//Functions to set properties of the system
	void InitializeFromStaticState(CStaticState<Dim> const &s_in, bool MakeUnstressed, bool UseUnitStiffness);
	void InitializeSimpleNetwork(vector<dvec> const &pos, vector< std::pair<int,int> > const &bonds, dbl _V);

//Set box
	void SetBox(CBox<Dim> *t_Box);				//!<Set the box
	void SetBoxPeriodic(int NonPeriodicDim=0);	//!<Set the box to be periodic
	void AssumeRectangularSymmetry();			//!<Set the assumption that the box is rectangular

	void SetPositions		(const Eigen::VectorXd &t_Positions);		//!<Set node positions
	void SetPositionsVirtual(const Eigen::VectorXd &t_Positions);
	void SetNodePosition		(const dvec &t_Position,int i);		//!<Set an individual node's position
	void SetNodePositionVirtual	(const dvec &t_Position,int i);
	
	void MoveNodes(const Eigen::VectorXd &t_Displacement);			//!<Move the nodes
	void MoveNodesVirtual(const Eigen::VectorXd &t_Displacement);

	void SetBondStiffness(dbl k0, int bi);
	void SetBondELength(dbl dr0, int bi);

	void SetVolume(dbl V);							//!<Set the volume without changing the radii

	void RemoveBonds(vector<bool> const &BondsToRemove);

//Functions to get properties of the system
	void GetPositions		(Eigen::VectorXd &) const;		//!<Copy the vector of positions
	void GetPositionsVirtual(Eigen::VectorXd &) const;
	void GetNodePosition		(dvec &, int i) const;	//!<Get an individual node's position
	void GetNodePositionVirtual	(dvec &, int i) const;
	void GetDisplacement		(int i, int j, dvec &displacement) const;	//!<Get the distence between two nodes
	void GetDisplacementVirtual	(int i, int j, dvec &displacement) const;

	void GetBondisBondjs(Eigen::VectorXi &, Eigen::VectorXi &) const;
	void Getij(int bi, int &i, int &j) const;
	int  Geti(int bi) const;
	int  Getj(int bi) const;

	void GetStiffnesses(Eigen::VectorXd &) const;		//!<Copy the vector of stiffnesses
	dbl  GetBondStiffness(int bi) const;
	
	void GetELengths(Eigen::VectorXd &) const;		//!<Copy the vector of ELengths
	dbl  GetBondELength(int bi) const;

	dbl  GetAverageELength() const;
	void GetBondInfo(int bi, int &i, int &j, dbl &dr0, dbl &k0) const;

	CBox<Dim>*  GetBox() const;			//!<Return a pointer to the box object
	CHarmonicSpringPotential* GetPotential() const;	//!<Return a pointer to the potential object
	
	dbl  GetVolume() const;				//!<Get the volume
	int  GetNodeNumber() const;
	int  GetNBonds() const;

	friend class CNetworkDatabase<Dim>;
};

template<int Dim>
CNetworkState<Dim>::CNetworkState()
	: N(0),
	  Nbonds(0),
	  Box(NULL),
	  Potential(NULL)
{}

/*
template<int Dim>
CNetworkState<DIM>::CNetworkState(int _N, int _Nb, dbl V)
	: N(_N),
	  Nbonds(_Nb),
	  Box(NULL),
	  Potential(NULL)
{
	SetPeriodicBox(Dim); //initialize to have free bcs.
	SetVolume(_V);

	if(Potential!=NULL)	delete Potential;
	Potential   = new CHarmonicSpringPotential(); //WARNING!! NOT GENERAL. This is the only potential I have implemented so far for springs

	Positions   = Eigen::VectorXd::Zero(Dim*N);
	Bondi       = Eigen::VectorXi::Zero(Nbonds);
	Bondj       = Eigen::VectorXi::Zero(Nbonds);
	Stiffnesses	= Eigen::VectorXd::Zero(Nbonds);
	ELengths    = Eigen::VectorXd::Zero(Nbonds);
}
*/

template<int Dim>
CNetworkState<Dim>::CNetworkState(CStaticState<Dim> const &s_in, bool MakeUnstressed, bool UseUnitStiffness)
	: N(0),
	  Nbonds(0),
	  Box(NULL),
	  Potential(NULL)
{
	InitializeFromStaticState(s_in, MakeUnstressed, UseUnitStiffness);
}

template<int Dim>
CNetworkState<Dim>::CNetworkState(CNetworkState<Dim> const &copy)
	: N(0),
	  Nbonds(0),
	  Box(NULL),
	  Potential(NULL)
{
	(*this) = copy;
}

template<int Dim>
static void cell_to_index(Eigen::Matrix<int,Dim,1> copies_per_side, Eigen::Matrix<int,Dim,1> cell, int &index)
{
	index = 0;
	int m = 1;
	for(int dd=0; dd<Dim; ++dd)
	{
		index += cell[dd]*m;
		m *= copies_per_side[dd];
	}
}
template<int Dim>
static void calculate_cell(int c, Eigen::Matrix<int,Dim,1> copies_per_side, Eigen::Matrix<int,Dim,1> &cell)
{
	cell = Eigen::Matrix<int,Dim,1>::Zero();
	for(int i=0; i<c; ++i)
	{
		++cell[0];
		for(int dd=0; dd<Dim; ++dd)
		{
			if(cell[dd] >= copies_per_side[dd])
			{
				assert( cell[dd]==copies_per_side[dd] );
				assert( dd < Dim-1 );
				cell[dd] -= copies_per_side[dd];
				++cell[dd+1];
			}
		}
	}
}

/*
template<int Dim>
static void calculate_offset(int c, Eigen::Matrix<int,Dim,1> copies_per_side, Eigen::Matrix<dbl,Dim,1> &offset)
{
	Eigen::Matrix<int,Dim,1> cell;
	calculate_cell(c, copies_per_side, cell);
	for(int dd=0; dd<Dim; ++dd)
		offset[dd] = ((dbl)cell[dd])/((dbl)copies_per_side[dd]);
}
*/

template<int Dim>
CNetworkState<Dim>::CNetworkState(CNetworkState<Dim> const &copy, ivec copies_per_side)
	: N(0),
	  Nbonds(0),
	  Box(NULL),
	  Potential(NULL)
{
	int Ncopies = 1;
	for(int dd=0; dd<Dim; ++dd)
	{
		assert(copies_per_side[dd]>0);
		Ncopies *= copies_per_side[dd];
	}

	if(Box!=NULL)		delete Box;
	if(Potential!=NULL)	delete Potential;
	
	Box			= copy.GetBox()->Clone();		//Clone() gives a deep copy
	Potential   = new CHarmonicSpringPotential(); //WARNING!! NOT GENERAL. This is the only potential I have implemented so far for springs

	//Adjust the size of the box.
	dmat Trans, Trans_new;
	Box->GetTransformation(Trans);
	dvec cps_dbl;
	for(int dd=0; dd<Dim; ++dd)	cps_dbl[dd] = (dbl)copies_per_side[dd];
	Trans_new = Trans * (cps_dbl.asDiagonal());
	Box->SetTransformation(Trans_new);
	
	///////////////////////////
	/////Set the positions/////
	///////////////////////////
	int Nsmall = copy.GetNodeNumber();
	N = Ncopies*Nsmall;
	Positions = Eigen::VectorXd::Zero(Dim*N);

	Eigen::VectorXd ptemp;
	copy.GetPositionsVirtual(ptemp);
	assert(ptemp.size()*Ncopies == Positions.size());

	//divide ptemp by copies_per_cell
	int ii, dd;
	for(ii=0; ii<Nsmall; ++ii)
		for(dd=0; dd<Dim; ++dd)
			ptemp[Dim*ii+dd] /= ((dbl)copies_per_side[dd]);

	dvec offset;
	int start;
	for(int c=0; c<Ncopies; ++c)
	{
		calculate_offset(c, copies_per_side, offset);
		start = Nsmall*c;
		for(ii=0; ii<Nsmall; ++ii)
		{
			for(dd=0; dd<Dim; ++dd)
				Positions[Dim*(start+ii) + dd] = ptemp[Dim*ii+dd]+offset[dd];
//			SetRadius(start+ii, copy.GetRadius(ii));
		}
	}

	///////////////////////
	/////Set the bonds/////
	///////////////////////
	int NBsmall = copy.GetNBonds();
	Nbonds = Ncopies*NBsmall;
	Bondi		= Eigen::VectorXi::Zero(Nbonds);
	Bondj		= Eigen::VectorXi::Zero(Nbonds);
	Stiffnesses	= Eigen::VectorXd::Zero(Nbonds);
	ELengths	= Eigen::VectorXd::Zero(Nbonds);

	int i,j;
	dbl k,l;
	dvec Displacement, DispTemp, diff;
	int jj;
	vector<bool> bondSet(Nbonds, false);
	for(int bi=0; bi<NBsmall; ++bi)
	{
		copy.GetBondInfo(bi, i, j, l, k);
		copy.GetDisplacement(i,j,Displacement);

		for(int c1=0; c1<Ncopies; ++c1)
		{
			ii = Nsmall*c1 + i;
			bool found = false;
			for(int c2=0; c2<Ncopies; ++c2)
			{
				jj = Nsmall*c2 + j;
				GetDisplacement(ii, jj, DispTemp);
				diff = Displacement - DispTemp;
				if(diff.norm() < 1e-8)
				{
					//add the bond
					int bbi = NBsmall*c1 + bi;
					assert( bondSet[bbi] == false );
					Bondi[bbi] = ii;
					Bondj[bbi] = jj;
					Stiffnesses[bbi] = k;
					ELengths[bbi] = l;
					bondSet[bbi] = true;
				
					found = true;
					break;
				}
			}
			if(!found) printf("WARNING! pair not found! %i %i\n", i, j);
		}
	}

/*

		//figure out if the bond crosses a boundary.
		ivec cell_offset;

		copy.GetBondInfo(bi, i, j, l, k);
		assert(i<Nsmall);
		assert(j<Nsmall);
		GetDisplacement(i,j,Displacement); //Displacement is the distance in the large system

		//Displacement is the distance from j to i (not the other way around!)
		for(int dd=0; dd<Dim; ++dd)
		{
			if(Displacement[dd] > 0.5/((dbl)copies_per_side[dd]))
				cell_offset[dd] = 1;
			else if(Displacement[dd] < -0.5/((dbl)copies_per_side[dd]))
				cell_offset[dd] = -1;
			else
				cell_offset[dd] = 0;
		}
		
		vector<bool> bondSet(Nbonds, false);
		for(int c=0; c<Ncopies; ++c)
		{
			ivec cell;
			calculate_cell(c, copies_per_side, cell);
			ivec nbr_cell = cell + cell_offset;

			//apply periodic boundary conditions
			for(int dd=0; dd<Dim; ++dd)
			{
				if(nbr_cell[dd] >= copies_per_side[dd])
					nbr_cell[dd] -= copies_per_side[dd];
				if(nbr_cell[dd] < 0)
					nbr_cell[dd] += copies_per_side[dd];
			}

			//Now connect node i in cell "cell" with node j in cell "nbr_cell"
			int nbr_cell_index;
			cell_to_index(copies_per_side, nbr_cell, nbr_cell_index);
			int ii, jj;
			ii = Nsmall*c + i;
			jj = Nsmall*nbr_cell_index + j;

			dvec disp1, disp2;
			GetDisplacement(ii,jj,disp1); //Displacement is the distance in the large system
			copy.GetDisplacement(i,j,disp2);
			dvec diff = disp1-disp2;
			if(diff.norm() > 1e-5)
			{
				printf("%e %e, i %i, j %i, %i %i\n", diff[0], diff[1], i, j, cell_offset[0], cell_offset[1]);
			}


			//add the bond
			int bbi = NBsmall*c + bi;
			assert( bondSet[bbi] == false );
			Bondi[bbi] = ii;
			Bondj[bbi] = jj;
			Stiffnesses[bbi] = k;
			ELengths[bbi] = l;
			bondSet[bbi] = true;
		}
	}
	*/
}	
	
	
template<int Dim>
void CNetworkState<Dim>::InitializeFromStaticState(CStaticState<Dim> const &s_in, bool MakeUnstressed, bool UseUnitStiffness)
{
	if(Box!=NULL)		delete Box;
	if(Potential!=NULL)	delete Potential;
	Box			= s_in.GetBox()->Clone();		//Clone() gives a deep copy
	Potential   = new CHarmonicSpringPotential(); //WARNING!! NOT GENERAL. This is the only potential I have implemented so far for springs

	CStaticState<Dim> s(s_in);
	CStaticComputer<Dim> c(s);
	if(c.StdPrepareSystem())
	{
		printf("Warning!!!!\n");
		assert(false);
	}
	c.CalculateStdData(false,false);

	N = c.Data.NPp;
	Nbonds = c.Bonds.GetNBonds();

	Positions   = Eigen::VectorXd::Zero(Dim*N);
	Bondi       = Eigen::VectorXi::Zero(Nbonds);
	Bondj       = Eigen::VectorXi::Zero(Nbonds);
	Stiffnesses	= Eigen::VectorXd::Zero(Nbonds);
	ELengths    = Eigen::VectorXd::Zero(Nbonds);
	assert(N == c.RattlerMap.size());

	//Set node positions
	Eigen::VectorXd PosTemp;
	s_in.GetPositionsVirtual(PosTemp);
	for(int im=0; im<N; ++im)
		for(int dd=0; dd<Dim; ++dd)
			Positions[Dim*im+dd] = PosTemp[Dim*c.RattlerMap[im]+dd];

	//Set equilibrium lengths and stiffnesses
	for(int bi=0; bi<Nbonds; ++bi)
	{
		Bondi[bi] = c.Bonds[bi].i;
		Bondj[bi] = c.Bonds[bi].j;

		dbl radsum = s.GetRadius(c.RattlerMap[c.Bonds[bi].i]) + s.GetRadius(c.RattlerMap[c.Bonds[bi].j]);
		if(MakeUnstressed)
			ELengths[bi] = c.Bonds[bi].rlen; //Use current distance between nodes
		else
			ELengths[bi] = radsum; //Use the sum of the radii of the original packing.

		if(UseUnitStiffness)
			Stiffnesses[bi] = 1.;
		else
			Stiffnesses[bi] = 1./POW2(radsum); //Not sure if this works for alpha!=2
	}
}

template<int Dim>
void CNetworkState<Dim>::InitializeSimpleNetwork(vector<dvec> const &pos, vector< std::pair<int,int> > const &bonds, dbl _V)
{
	SetBoxPeriodic(Dim); //initialize to have free bcs.
	SetVolume(_V);

	if(Potential!=NULL)	delete Potential;
	Potential   = new CHarmonicSpringPotential(); //WARNING!! NOT GENERAL. This is the only potential I have implemented so far for springs

	N = (int)pos.size();
	Nbonds = bonds.size();

	Positions   = Eigen::VectorXd::Zero(Dim*N);
	Bondi       = Eigen::VectorXi::Zero(Nbonds);
	Bondj       = Eigen::VectorXi::Zero(Nbonds);
	Stiffnesses	= Eigen::VectorXd::Zero(Nbonds);
	ELengths    = Eigen::VectorXd::Zero(Nbonds);

	for(int i=0; i<N; ++i)
		SetNodePosition(pos[i], i);
//		for(int dd=0; dd<Dim; ++dd)
//			Positions[Dim*i+dd] = pos[i][dd];

	int ii, jj;
	dvec disp;
	for(int bi=0; bi<Nbonds; ++bi)
	{
		ii = bonds[bi].first;
		jj = bonds[bi].second;
		assert(ii >= 0);
		assert(jj >= 0);
		assert(ii < N);
		assert(jj < N);
		Bondi[bi] = ii;
		Bondj[bi] = jj;

		//make bond unstressed
		GetDisplacement(ii,jj,disp);
		ELengths[bi] = disp.norm();

		//set unit stiffness
		Stiffnesses[bi] = 1.;
	}

}





template<int Dim>	
CNetworkState<Dim>::~CNetworkState()
{
	delete Box;
	delete Potential;
}

template<int Dim>
CNetworkState<Dim> &CNetworkState<Dim>::operator=(CNetworkState<Dim> const &copy)
{
	if(this != &copy)
	{
		N = copy.GetNodeNumber();
		Nbonds = copy.GetNBonds();
	
		if(Box!=NULL)		delete Box;
		if(Potential!=NULL)	delete Potential;
		
		Box			= copy.GetBox()->Clone();		//Clone() gives a deep copy
		Potential	= copy.GetPotential()->Clone();	//Clone() gives a deep copy
		
		copy.GetPositionsVirtual(Positions);
		copy.GetBondisBondjs(Bondi,Bondj);
		copy.GetStiffnesses(Stiffnesses);
		copy.GetELengths(ELengths);
	}
	return *this;
}



template <int Dim>
void CNetworkState<Dim>::SetBox(CBox<Dim> *t_Box)
{
	if(Box!=NULL)	delete Box;
	Box = t_Box->Clone();
}

template<int Dim>
void CNetworkState<Dim>::SetBoxPeriodic(int NonPeriodicDim)
{
	if(Box!=NULL)	delete Box;
	switch(NonPeriodicDim)
	{
		case 0:		Box = new CPeriodicBox<Dim,0>();	break;
		case 1:		Box = new CPeriodicBox<Dim,1>();	break;
		case Dim:	Box = new CPeriodicBox<Dim,Dim>();	break;
		default:	assert(false);
	}
}

template<int Dim>
void CNetworkState<Dim>::AssumeRectangularSymmetry()
{
	Box->SetSymmetry(Box->RECTANGULAR);
}






template <int Dim>
void CNetworkState<Dim>::SetPositions(const Eigen::VectorXd &t_Positions)
{
	//set the position vector
	if(t_Positions.rows()!=Dim*N)
		throw(CException("CNetworkState<Dim>::SetPositions","Positions vector has an inconsistent size.")); 
	
	Positions = t_Positions;

	Box->InverseTransform(Positions);
}

template <int Dim>
void CNetworkState<Dim>::SetPositionsVirtual(const Eigen::VectorXd &t_Positions)
{
	//set the position vector
	if(t_Positions.rows()!=Dim*N)
		throw(CException("CNetworkState<Dim>::SetPositionsVirtual","Positions vector has an inconsistent size.")); 
		
	Positions = t_Positions;
}

template <int Dim>
void CNetworkState<Dim>::SetNodePosition(const Eigen::Matrix<dbl,Dim,1> &t_Positions,int i)
{
	if(i<0||i>N)
		throw(CException("CNetworkState<Dim>::SetNodePosition","Attempting to set the position of a particle that does not exist."));
	
	Eigen::Matrix<dbl,Dim,1> tpos = t_Positions;
	Box->InverseTransform(tpos);
	Positions.segment<Dim>(Dim*i) = tpos;
}

template <int Dim>
void CNetworkState<Dim>::SetNodePositionVirtual(const Eigen::Matrix<dbl,Dim,1> &t_Positions,int i)
{
	if(i<0||i>N)
		throw(CException("CNetworkState<Dim>::SetNodePositionVirtual","Attempting to set the position of a particle that does not exist."));
		
	Positions.segment<Dim>(Dim*i) = t_Positions;
}

template<int Dim>
void CNetworkState<Dim>::MoveNodes(const Eigen::VectorXd &t_Displacement)
{
	if(t_Displacement.rows()!=Dim*N)
		throw(CException("CNetworkState<Dim>::MoveNodes","Displacement vector has an inconsistent size.")); 
		
//	Box->InverseTransform(t_Displacement);
//	Box->MoveNodes(Positions,t_Displacement);
	Box->InverseTransformAndMove(Positions,t_Displacement);
}

template<int Dim>
void CNetworkState<Dim>::MoveNodesVirtual(const Eigen::VectorXd &t_Displacement)
{
	if(t_Displacement.rows()!=Dim*N)
		throw(CException("CNetworkState<Dim>::MoveNodesVirtual","Displacement vector has an inconsistent size.")); 
		
	Box->MoveParticles(Positions,t_Displacement);
}

template<int Dim>
void CNetworkState<Dim>::SetBondStiffness(dbl k0, int bi)
{
	assert(bi>=0&&bi<Nbonds);
	Stiffnesses[bi] = k0;
}

template<int Dim>
void CNetworkState<Dim>::SetBondELength(dbl dr0, int bi)
{
	assert(bi>=0&&bi<Nbonds);
	ELengths[bi] = dr0;
}


template <int Dim>
void CNetworkState<Dim>::SetVolume(dbl V)
{
	Box->SetVolume(V);
}


/**
 *  Bonds that are not removed are not altered in any way, and N does not change, nor does Positions.
 *
 *  @param[in] BondsToRemove A constant vector<bool> of length list.size() that determines
 *  which bonds to remove. Bond i is removed only if BondsToRemove[i]==true.
 */
template<int Dim>
void CNetworkState<Dim>::RemoveBonds(vector<bool> const &BondsToRemove)
{
	assert(BondsToRemove.size() == Nbonds);
	assert(Nbonds == Bondi.size());
	assert(Nbonds == Bondj.size());
	assert(Nbonds == Stiffnesses.size());
	assert(Nbonds == Stiffnesses.size());

	//First, count the number of bonds that will be left.
	int numBondsLeft = 0;
	for(vector<bool>::const_iterator it=BondsToRemove.begin(); it!=BondsToRemove.end(); ++it)
		if( (*it) == 0 ) ++numBondsLeft;
	assert(numBondsLeft <= Nbonds);
	assert(numBondsLeft >= 0);
	if(numBondsLeft == Nbonds) return;

	//Now, make new copies of Bondi, Bondj, Stiffnesses, and ELengths that omit the removed bonds
	Eigen::VectorXi t_Bondi			= Eigen::VectorXi::Zero(numBondsLeft);
	Eigen::VectorXi t_Bondj			= Eigen::VectorXi::Zero(numBondsLeft);
	Eigen::VectorXd t_Stiffnesses	= Eigen::VectorXd::Zero(numBondsLeft);
	Eigen::VectorXd t_ELengths		= Eigen::VectorXd::Zero(numBondsLeft);

	int index = 0;
	for(int ii=0; ii<Nbonds; ++ii)
	{
		if( BondsToRemove[ii] == 0 )
		{
			t_Bondi[index]			= Bondi[ii];
			t_Bondj[index]			= Bondj[ii];
			t_Stiffnesses[index]	= Stiffnesses[ii];
			t_ELengths[index]		= ELengths[ii];

			++index;
		}
	}
	assert(index == numBondsLeft);

	//Finally, replace the old vectors with the new vectors, and change Nbonds
	Bondi		= t_Bondi;
	Bondj		= t_Bondj;
	Stiffnesses	= t_Stiffnesses;
	ELengths	= t_ELengths;
	Nbonds = numBondsLeft;

	assert(Nbonds == Bondi.size());
	assert(Nbonds == Bondj.size());
	assert(Nbonds == Stiffnesses.size());
	assert(Nbonds == Stiffnesses.size());

/*
	assert(BondsToRemove.size() == list.size());
	vector<BOND> tlist;
	tlist.reserve(list.size());

	for(int i=0; i<(int)list.size(); ++i)
	if(!BondsToRemove[i])
	tlist.push_back(list[i]);

	tlist.swap(list);
	*/
}






template<int Dim>
void CNetworkState<Dim>::GetPositions(Eigen::VectorXd &t_Positions) const
{
	t_Positions = Positions;
	Box->Transform(t_Positions);
}

template<int Dim>
void CNetworkState<Dim>::GetPositionsVirtual(Eigen::VectorXd &t_Positions) const
{
	t_Positions = Positions;
}

template <int Dim>
void CNetworkState<Dim>::GetNodePosition(Eigen::Matrix<dbl,Dim,1> &t_Position, int i) const
{
	if(i<0||i>N)
		throw(CException("CNetworkState<Dim>::GetNodePosition","Attempting to set the position of a particle that does not exist."));

	t_Position = Positions.segment<Dim>(Dim*i);
	Box->Transform(t_Position);
}

template <int Dim>
void CNetworkState<Dim>::GetNodePositionVirtual(Eigen::Matrix<dbl,Dim,1> &t_Position, int i) const
{
	if(i<0||i>N)
		throw(CException("CNetworkState<Dim>::GetNodePositionVirtual","Attempting to set the position of a particle that does not exist."));

	t_Position = Positions.segment<Dim>(Dim*i);
}

template <int Dim>
void CNetworkState<Dim>::GetDisplacement(int i,int j,Eigen::Matrix<dbl,Dim,1> &displacement) const
{
	if(i<0||i>N||j<0||j>N)
		throw(CException("CNetworkState<Dim>::GetDisplacement","Attempting to get the displacement between particles that do not exist."));

	Box->MinimumDisplacement(Positions.segment<Dim>(Dim*i),Positions.segment<Dim>(Dim*j),displacement);
	
	Box->Transform(displacement);
}

template <int Dim>
void CNetworkState<Dim>::GetDisplacementVirtual(int i,int j,Eigen::Matrix<dbl,Dim,1> &displacement) const
{
	if(i<0||i>N||j<0||j>N)
		throw(CException("CNetworkState<Dim>::GetDisplacement","Attempting to get the displacement between particles that do not exist."));

	Box->MinimumDisplacement(Positions.segment<Dim>(Dim*i),Positions.segment<Dim>(Dim*j),displacement);
}


template <int Dim>
void CNetworkState<Dim>::GetBondisBondjs(Eigen::VectorXi &t_Bondi, Eigen::VectorXi &t_Bondj) const
{
	t_Bondi = Bondi;
	t_Bondj = Bondj;
}

template <int Dim>
void CNetworkState<Dim>::Getij(int bi, int &i, int &j) const
{
	i = Bondi[bi];
	j = Bondj[bi];
}

template <int Dim>
int CNetworkState<Dim>::Geti(int bi) const
{
	return Bondi[bi];
}

template <int Dim>
int CNetworkState<Dim>::Getj(int bi) const
{
	return Bondj[bi];
}

template <int Dim>
void CNetworkState<Dim>::GetStiffnesses(Eigen::VectorXd &t_Stiffnesses) const
{
	t_Stiffnesses = Stiffnesses;
}

template <int Dim>
dbl CNetworkState<Dim>::GetBondStiffness(int bi) const
{
	return Stiffnesses[bi];
}

template <int Dim>
void CNetworkState<Dim>::GetELengths(Eigen::VectorXd &t_ELengths) const
{
	t_ELengths = ELengths;
}

template <int Dim>
dbl CNetworkState<Dim>::GetBondELength(int bi) const
{
	return ELengths[bi];
}


template <int Dim>
dbl CNetworkState<Dim>::GetAverageELength() const
{
	return ELengths.mean();
}

template <int Dim>
void CNetworkState<Dim>::GetBondInfo(int bi, int &i, int &j, dbl &dr0, dbl &k0) const
{
	i   = Bondi[bi];
	j   = Bondj[bi];
	dr0 = ELengths[bi];
	k0  = Stiffnesses[bi];
}

template <int Dim>
CBox<Dim> *CNetworkState<Dim>::GetBox() const
{
	return Box;
}

template <int Dim>
CHarmonicSpringPotential *CNetworkState<Dim>::GetPotential() const
{
	return Potential;
}


template <int Dim>
dbl CNetworkState<Dim>::GetVolume() const
{
	return Box->CalculateVolume();
}
	
template <int Dim>
int CNetworkState<Dim>::GetNodeNumber() const
{
	return N;
}

template <int Dim>
int CNetworkState<Dim>::GetNBonds() const
{
	return Nbonds;
}














}


#endif //NETWORK_STATE
