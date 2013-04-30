#ifndef STATIC_STATE

#define STATIC_STATE

/////////////////////////////////////////////////////////////////////////////////
//Static state class. 
//
//Description
//		This class stores information to represent systems of isotropic particles
//		and implements functions to manipulate this information. The data that is
//		fundamentally necessary to represent this system are particle positions
//		and radii which are stored as Eigen vectors. Other information that is 
//		required is a potential energy function and boundary conditions. These 
//		ideas are implemented through inherited classes. Writing to and reading from 
//		NetCDF files is implemented.
//
//		Particle positions are stored in the unit box [0,1]x...x[0,1] and then mapped
//		into a system of proper shape by the box class. The stored positions are
//		denoted "virtual" positions whereas the actual positions are denoted "real".
//
//Variables
// 		Particle positions as an eigen vector
//		Particle radii as an eigen vector
//		Number of particles as an integer
//		Box as a CBox class
//		Potential as a CPotential class
//		
//Implements
//		Reading/Writing to netcdf files.
//		Getting/Setting particle positions.
//						particle radii.
//						box.
//						potential.
//
//
//File Formats
//		NetCDF
//			## Note: only files of a single N and Dimension may be stored in a given 
//			## file. To denote whether the NetCDF file has been populated with variables,
//			## dimensions, etc... an attribute "System_Populated" gets created.
//			-> Three dimensions: record number (unlimited), degrees of freedom (DOF), dimension,
//			   and particle number.
//			-> Particle positions are stored as a variable.
//			-> Particle radii are stored as a variable.
//			-> Box is stored as a box (see CBox definition.)
//			-> Potential is stored as a potential (see CPotential definition.)
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include "../Potentials/Potential.h"
#include "../Boundaries/Box.h"
#include "../Potentials/HarmonicPotential.h"
#include "../Boundaries/PeriodicBox.h"
#include "netcdfcpp.h"
#include "../Resources/MersenneTwister.h"

namespace LiuJamming
{

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
class CStaticState
{
private:
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;

	Eigen::VectorXd Positions;
	Eigen::VectorXd Radii;

	int N;

	CBox<Dim> *Box;
	CPotential *Potential;
		
//Functions to check that the correct dimensions, variables, and attributes are present
//in a netCDF file.
	bool CheckNetCDF(const NcFile &file);
	void PopulateNetCDF(NcFile &file);
	
public:
//Constructors/Destructors and copy operators
	CStaticState();
	CStaticState(int _N);
	CStaticState(const CStaticState &copy);
	
	~CStaticState();

	CStaticState<Dim> &operator=(const CStaticState<Dim> &copy);

//Functions to construct systems 
	void RandomizePositions(long seed = 1);
	void RandomizePositions(MTRand *random);
	void Read(const NcFile &File,int Record);
	
//Functions to write systems
	void Write(NcFile &File,int Record);

//Functions to set properties of the system
	void SetRadii(const Eigen::VectorXd &t_Radii);
	void SetRadius(int i,dbl r);
	void SetPositions(const Eigen::VectorXd &t_Positions);
	void SetPositionsVirtual(const Eigen::VectorXd &t_Positions);
	void SetParticlePosition(const dvec &t_Positions,int i);
	void SetParticlePositionVirtual(const dvec &t_Positions,int i);
	void MoveParticles(const Eigen::VectorXd &t_Displacement);	
	void MoveParticlesVirtual(const Eigen::VectorXd &t_Displacement);
	void SetBox(CBox<Dim> *t_Box);
	void SetPotential(CPotential *t_Potential);
	void SetPackingFraction(dbl phi);

//Functions to get properties of the system
	void GetRadii(Eigen::VectorXd &);	
	dbl GetRadius(int i);
	void GetPositions(Eigen::VectorXd &);
	void GetPositionsVirtual(Eigen::VectorXd &);
	void GetParticlePosition(dvec &, int i);
	void GetParticlePositionVirtual(dvec &, int i);
	void GetDisplacement(int i, int j, dvec &displacement);
	void GetDisplacementVirtual(int i, int j, dvec &displacement);
	CBox<Dim> *GetBox();	
	CPotential *GetPotential();
	dbl  GetSphereVolume() const;
	dbl  GetPackingFraction() const;
	dbl  GetVolume() const;
	int  GetParticleNumber() const;

	dvec GetMaxDistance() const;

	void PrintParticles() const;

//	//Allow the class CGrid<Dim> to access private information.
//	friend class CGrid<Dim>;
};


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


//Functions to check that the correct dimensions, variables, and attributes are present
//in a netCDF file.
template <int Dim>
bool CStaticState<Dim>::CheckNetCDF(const NcFile &file)
{
	return (file.get_att("System_Populated")!=NULL);
}

template <int Dim>
void CStaticState<Dim>::PopulateNetCDF(NcFile &file)
{
	NcDim *records = file.get_dim("Records");
	if(records==NULL);
		records = file.add_dim("Records");
		
	NcDim *DOF = file.add_dim("System_DOF",Dim*N);
	NcDim *Number = file.add_dim("System_Number",N);
	file.add_dim("System_Dimension",Dim);
	
	file.add_att("System_Populated",1);
	
	file.add_var("System_Positions",ncDouble,records, DOF);
	file.add_var("System_Radii",ncDouble,records, Number);
}

//Constructors/Destructors and copy operators
//
//!!!POTENTIAL MEMORY LEAK?
template <int Dim>
CStaticState<Dim>::CStaticState()
{
	N = 0;
	Box = new CPeriodicBox<Dim>();
	Potential = new CHarmonicPotential();
}

template <int Dim>
CStaticState<Dim>::CStaticState(int _N)
{
	N = _N;
	Box = new CPeriodicBox<Dim>();
	Potential = new CHarmonicPotential();
	
	Positions = Eigen::ArrayXd::Zero(Dim*N);
	Radii = Eigen::ArrayXd::Constant(N,1.);
//	Radii = Eigen::ArrayXd::Constant(N,0.5);
}
	
template<int Dim>
CStaticState<Dim>::CStaticState(const CStaticState &copy)
{
	N = copy.GetParticleNumber();
	
	Box = copy.GetBox()->copy(); //????????
	Potential = copy.GetPotential()->copy();
	
	Positions = copy.GetPositions();
	Radii = copy.GetRadii();
}

template<int Dim>	
CStaticState<Dim>::~CStaticState()
{
	delete Box;
	delete Potential;
}

template<int Dim>
CStaticState<Dim> &CStaticState<Dim>::operator=(const CStaticState<Dim> &copy)
{
	if(this != &copy)
	{
		N = copy.GetParticleNumber();
		
		Box = copy.GetBox()->copy();
		Potential = copy.GetPotential()->copy();
		
		Positions = copy.GetPositions();
		Radii = copy.GetRadii();
	}
	return *this;
}

//Functions to construct systems 
//Randomly positions particles in the system. 
template<int Dim>
void CStaticState<Dim>::RandomizePositions(MTRand *random)
{
	for(int i = 0; i<N ; i++)
		for(int d = 0 ; d < Dim; d++)
			Positions(Dim*i + d) = random->rand();
}

template<int Dim>
void CStaticState<Dim>::RandomizePositions(long seed)
{
	MTRand *random = new MTRand(seed);
	RandomizePositions(random);
	delete random;
}

//Read a system from a netcdf file.
template <int Dim>
void CStaticState<Dim>::Read(const NcFile &File,int Record)
{
	//Check to make sure the netcdf file is populated
	if(!CheckNetCDF(File))
		throw(CException("CStaticState<Dim>::ReadSystem","Attempting to read a system from a file that has no system information.")); 
	
	//check to make sure the number of particles and dimension are consistent
	if(File.get_dim("System_Number")->size()!=N||File.get_dim("System_Dimension")->size()!=Dim)
		throw(CException("CStaticState<Dim>::ReadSystem","Attempting to read a system from a file that has inconsistent system information.")); 

	//check to make sure the record number is within the bounds
	if(Record>=File.get_dim("Records")->size())
		throw(CException("CStatciState<Dim>::ReadSystem","Attempting to read a record that does not exist.")); 

	//read the position and radius information
	NcVar *PositionsVar = File.get_var("System_Positions");		
	NcVar *RadiiVar = File.get_var("System_Radii");
	
	PositionsVar->set_cur(Record);
	PositionsVar->get(Positions.data(),1,Dim*N);
	
	RadiiVar->set_cur(Record);
	RadiiVar->get(Radii.data(),1,N);
	
	//read the box and potential information.
	if(Box!=NULL)
		delete Box;
	Box = CBox<Dim>::Read(File,Record);
	
	if(Potential!=NULL)
		delete Potential;
	Potential = CPotential::Read(File,Record);

}
	
//Functions to write systems
template <int Dim>
void CStaticState<Dim>::Write(NcFile &File,int Record)
{
	//If the netcdf file is not populated then populate it
	if(!CheckNetCDF(File))
		PopulateNetCDF(File);
	
		
	//check to make sure the number of particles and dimension are consistent
	if(File.get_dim("System_Number")->size()!=N||File.get_dim("System_Dimension")->size()!=Dim)
		throw(CException("CStaticState<Dim>::Write","Attempting to write a system from a file that has inconsistent system information.")); 

	//check to make sure the record number is within the bounds
	if(Record>File.get_dim("Records")->size())
		throw(CException("CStaticState<Dim>::Write","Attempting to write a record that does not exist.")); 
		
	//write the position and radius information
	NcVar *PositionsVar = File.get_var("System_Positions");		
	NcVar *RadiiVar = File.get_var("System_Radii");
	
	cout << "Saving positions.\n";
	
	PositionsVar->set_cur(Record);
	PositionsVar->put(Positions.data(),1,Dim*N);
	
	cout << "Saving radii.\n";
	
	RadiiVar->set_cur(Record);
	RadiiVar->put(Radii.data(),1,N);
	
	cout << "Saving box.\n";
	//read the box and potential information.
	Box->Write(File,Record);
	cout << "Saving potential.\n";
	Potential->Write(File,Record);	
}
		
		
//Functions to set properties of the system
template <int Dim>
void CStaticState<Dim>::SetRadii(const Eigen::VectorXd &t_Radii)
{
	//set the radii vector  
	if(t_Radii.rows()!=N)
		throw(CException("CStaticState<Dim>::SetRadii","Radii vector has an inconsistent size.")); 

	Radii = t_Radii;
}

template <int Dim>
void CStaticState<Dim>::SetRadius(int i, dbl r)
{
	//set the radii vector  
	if(i>N)
		throw(CException("CStaticState<Dim>::SetRadius","Attempting to set the radius of a particle that doesn't exist.")); 

	Radii(i) = r;
}

template <int Dim>
void CStaticState<Dim>::SetPositions(const Eigen::VectorXd &t_Positions)
{
	//set the position vector
	if(t_Positions.rows()!=Dim*N)
		throw(CException("CStaticState<Dim>::SetPositions","Positions vector has an inconsistent size.")); 
	
	Positions = t_Positions;

	Box->InverseTransform(Positions);
}

template <int Dim>
void CStaticState<Dim>::SetPositionsVirtual(const Eigen::VectorXd &t_Positions)
{
	//set the position vector
	if(t_Positions.rows()!=Dim*N)
		throw(CException("CStaticState<Dim>::SetPositionsVirtual","Positions vector has an inconsistent size.")); 
		
	Positions = t_Positions;
}

template <int Dim>
void CStaticState<Dim>::SetParticlePosition(const Eigen::Matrix<dbl,Dim,1> &t_Positions,int i)
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::SetParticlePosition","Attempting to set the position of a particle that does not exist."));
		
	Box->InverseTransform(t_Positions);
	Positions.segment<Dim>(Dim*i) = t_Positions;
}

template <int Dim>
void CStaticState<Dim>::SetParticlePositionVirtual(const Eigen::Matrix<dbl,Dim,1> &t_Positions,int i)
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::SetParticlePositionVirtual","Attempting to set the position of a particle that does not exist."));
		
	Positions.segment<Dim>(Dim*i) = t_Positions;
}

template<int Dim>
void CStaticState<Dim>::MoveParticles(const Eigen::VectorXd &t_Displacement)
{
	if(t_Displacement.rows()!=Dim*N)
		throw(CException("CStaticState<Dim>::MoveParticles","Displacement vector has an inconsistent size.")); 
		
//	Box->InverseTransform(t_Displacement);
//	Box->MoveParticles(Positions,t_Displacement);
	Box->InverseTransformAndMove(Positions,t_Displacement);
}

template<int Dim>
void CStaticState<Dim>::MoveParticlesVirtual(const Eigen::VectorXd &t_Displacement)
{
	if(t_Displacement.rows()!=Dim*N)
		throw(CException("CStaticState<Dim>::MoveParticlesVirtual","Displacement vector has an inconsistent size.")); 
		
	Box->MoveParticles(Positions,t_Displacement);
}

template <int Dim>
void CStaticState<Dim>::SetBox(CBox<Dim> *t_Box)
{
	Box = t_Box;
}

template <int Dim>
void CStaticState<Dim>::SetPotential(CPotential *t_Potential)
{
	Potential = t_Potential;
}

template <int Dim>
void CStaticState<Dim>::SetPackingFraction(dbl phi)
{
	Box->SetVolume(GetSphereVolume()/phi);
	assert( fabs(GetPackingFraction()-phi) < 1e-10 );
}



//Functions to get properties of the system
template <int Dim>
void CStaticState<Dim>::GetRadii(Eigen::VectorXd &t_Radii)
{
	t_Radii = Radii;
}

template <int Dim>
dbl CStaticState<Dim>::GetRadius(int i)
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::GetRadius","Attempting to get the radius of a particle that does not exist."));
	return Radii(i);
}

template<int Dim>
void CStaticState<Dim>::GetPositions(Eigen::VectorXd &t_Positions)
{
	t_Positions = Positions;
	Box->Transform(t_Positions);
}

template<int Dim>
void CStaticState<Dim>::GetPositionsVirtual(Eigen::VectorXd &t_Positions)
{
	t_Positions = Positions;
}

template <int Dim>
void CStaticState<Dim>::GetParticlePosition(Eigen::Matrix<dbl,Dim,1> &t_Position, int i)
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::GetParticlePosition","Attempting to set the position of a particle that does not exist."));

	t_Position = Positions.segment<Dim>(Dim*i);
	Box->Transform(t_Position);
}

template <int Dim>
void CStaticState<Dim>::GetParticlePositionVirtual(Eigen::Matrix<dbl,Dim,1> &t_Position, int i)
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::GetParticlePositionVirtual","Attempting to set the position of a particle that does not exist."));

	t_Position = Positions.segment<Dim>(Dim*i);
}

template <int Dim>
void CStaticState<Dim>::GetDisplacement(int i,int j,Eigen::Matrix<dbl,Dim,1> &displacement)
{
	if(i<0||i>N||j<0||j>N)
		throw(CException("CStaticState<Dim>::GetDisplacement","Attempting to get the displacement between particles that do not exist."));

	Box->MinimumDisplacement(Positions.segment<Dim>(Dim*i),Positions.segment<Dim>(Dim*j),displacement);
	
	Box->Transform(displacement);
}

template <int Dim>
void CStaticState<Dim>::GetDisplacementVirtual(int i,int j,Eigen::Matrix<dbl,Dim,1> &displacement)
{
	if(i<0||i>N||j<0||j>N)
		throw(CException("CStaticState<Dim>::GetDisplacement","Attempting to get the displacement between particles that do not exist."));

	Box->MinimumDisplacement(Positions.segment<Dim>(Dim*i),Positions.segment<Dim>(Dim*j),displacement);
}

template <int Dim>
CBox<Dim> *CStaticState<Dim>::GetBox()
{
	return Box;
}

template <int Dim>
CPotential *CStaticState<Dim>::GetPotential()
{
	return Potential;
}

template <int Dim>
dbl CStaticState<Dim>::GetSphereVolume() const
{
	dbl SphereVolume = 0.;
	for(int i=0; i<Radii.size(); ++i)
		SphereVolume += std::pow(Radii[i],Dim);
	SphereVolume *= nSphere_Vn(Dim);
	return SphereVolume;
}

template <int Dim>
dbl CStaticState<Dim>::GetPackingFraction() const
{
	return GetSphereVolume()/GetVolume();
}

template <int Dim>
dbl CStaticState<Dim>::GetVolume() const
{
	return Box->CalculateVolume();
}
	
template <int Dim>
int CStaticState<Dim>::GetParticleNumber() const
{
	return N;
}

template <int Dim>
dvec CStaticState<Dim>::GetMaxDistance() const
{

	dbl MaxDistance = 2.*Radii.maxCoeff()*Potential->ComputeSupport();
	dvec = dvec::Constant(MaxDistance);

}
	
template <int Dim>
void CStaticState<Dim>::PrintParticles() const
{
	Eigen::VectorXd RealPos = Positions;
	Box->Transform(RealPos);
	for(int i=0; i<N; ++i)
	{
		printf("sphere:% 5i   Position: ", i); 
		for(int dd=0; dd<Dim; dd++) printf("% 20.14f  ", RealPos[Dim*i+dd]);
		printf("  Radius: %16.14f", Radii[i]);
		printf("\n");
	}
}

}

#endif
