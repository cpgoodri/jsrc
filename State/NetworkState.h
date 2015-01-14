#ifndef NETWORK_STATE
#define NETWORK_STATE

#include "../State/StaticState.h"
#include "../Computers/StaticComputer.h"
#include "../Computers/BondList.h"
#include "../Potentials/HarmonicSpringPotential.h"



namespace LiuJamming
{

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
	CNetworkState(CStaticState<Dim> const &s_in, bool MakeUnstressed, bool UseUnitStiffness);
	CNetowrkState(CNetworkState const &copy);								//!<Copy constructor
	CNetworkState<Dim> &operator=(CNetworkState<Dim> const &copy);			//!<Copy operator
	~CNetworkState();

//Functions to set properties of the system
	void InitializeFromStaticState(CStaticState<Dim> const &s_in, bool MakeUnstressed, bool UseUnitStiffness);

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
};

template<int Dim>
CNetworkState<Dim>::CNetworkState()
	: N(0),
	  Nbonds(0),
	  Box(NULL),
	  Potential(NULL)
{}

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
CNetowrkState<Dim>::CNetowrkState(CNetworkState<Dim> const &copy)
	: N(0),
	  Nbonds(0),
	  Box(NULL),
	  Potential(NULL)
{
	(*this) = copy;
}

/*
template<int Dim>
CNetworkState<Dim>::CNetworkState(CStaticState<Dim> const &s_in, bool MakeUnstressed, bool UseUnitStiffness)
	: N(0),
	  Box(NULL),
	  Potential(NULL)
{
	if(Box!=NULL)		delete Box;
	if(Potential!=NULL)	delete Potential;
	Box			= s_in.GetBox()->Clone();		//Clone() gives a deep copy
	Potential   = new CHarmonicSpringPotential(); //This is the only potential I have implemented so far for springs

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
*/

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
void CNetworkState<Dim>::Geti(int bi) const
{
	return Bondi[bi];
}

template <int Dim>
void CNetworkState<Dim>::Getj(int bi) const
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
