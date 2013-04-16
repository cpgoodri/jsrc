#ifndef GRID

#define GRID

/////////////////////////////////////////////////////////////////////////////////
//Grid class. 
//
//Description
//		This class maintains a spatial partition of particle systems so that
//		only local regions of a given particle need to be searched for neighbors.
//		The f 
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
//			-> Three dimensions: record number (unlimited), degrees of freedom (DOF), and
//			   particle number.
//			-> Dimension is stored as an attribute
//			-> Particle positions are stored as a variable.
//			-> Particle radii are stored as a variable.
//			-> Box is stored as a box (see CBox definition.)
//			-> Potential is stored as a potential (see CPotential definition.)
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include "../Potentials/Potential.h"
#include "../Boundaries/Box.h"
#include "../Resources/MersenneTwister.h"
#include <list>


namespace LiuJamming
{

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
class CGrid
{
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
public:

	class iterator {
	private:
		list<int>::iterator CurrentCell;
		list<int>::iterator End;
		
		int CurrentParticle;
		
		int *OccupancyList;
		int *CellList;
		
	public:
		iterator(list<int> &dual,int *occ_list, int *cell_list)
		{
			CurrentCell = dual.begin();
			End = dual.end();
			OccupancyList = occ_list;
			CellList = cell_list;
			CurrentParticle = CellList[(*CurrentCell)];
		}
	
		iterator &operator ++() {
			if(CurrentParticle!=-1){
				CurrentParticle = OccupancyList[CurrentParticle];
				if(CurrentParticle==-1&&CurrentCell!=End)
				{
					CurrentCell++;
					CurrentParticle = CellList[(*CurrentCell)];
				}
			}
			return (*this);
		}

		void operator ++(int dummy) {
			if(CurrentParticle!=-1){
				//cout << "Incrementing iterator; current particle is " << CurrentParticle << endl;
				CurrentParticle = OccupancyList[CurrentParticle];
				while(CurrentParticle==-1&&CurrentCell!=End)
				{
					CurrentCell++;
					//cout << "Incrementing cell iterator; current cell is " << (*CurrentCell) << endl;
					if(CurrentCell!=End)
						CurrentParticle = CellList[(*CurrentCell)];
				}
			}
		}
		
		int operator *()
		{
			return CurrentParticle;
		}
	};

private:
	//Box and particle information
	CStaticState<Dim> *State;
	int N;

	//Cell size given by twice the maximum radius.
	dvec CellSize;
	
	//The number of cells in each dimension
	int N_Cells[Dim];
	int TotalCells;
	
	//The list of particles in each cell
	int *OccupancyList;
	int *CellList;
	
	//The list of connections between cells.
	list<int> *DualList;
	
	//A boolean to indicate whether we have to reallocate the grid.
	bool Reallocate;
	
public:
//Constructor and copy operators
	CGrid(CStaticState<Dim> *_s);
	CGrid(const CGrid &copy);
	
	const CGrid<Dim> &operator=(const CGrid<Dim> &copy);
	
	void SetState(CStaticState<Dim> *s);

//Functions to construct the grid
	void Allocate();
	void Construct();

//Functions to access the grid
 	void CellToCoordinates(int i, dvec &coordinates);
	int CoordinatesToCell(const dvec &coordinates);
	iterator GetParticleIterator(int i);
	
//Function to compute the grid cell size ; 
	dvec ComputeCellSize();
	

};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//Constructor and copy operators
template <int Dim>
CGrid<Dim>::CGrid(CStaticState<Dim> *s) : State(s)
{
	N = State->GetParticleNumber();
	OccupancyList = new int[N];
	for(int i = 0 ; i<N ; i++)
		OccupancyList[i] = -1;
	CellList = NULL;
	DualList = NULL;
}

template <int Dim>
CGrid<Dim>::CGrid(const CGrid &copy) : State(copy.State)
{
	N = State->GetParticleNumber();
	OccupancyList = new int[N];
	for(int i = 0 ; i<N ; i++)
		OccupancyList[i] = -1;
	CellList = NULL;
	DualList = NULL;
}

//Possible memory leak!
template <int Dim>
const CGrid<Dim> &CGrid<Dim>::operator=(const CGrid<Dim> &copy)
{
	State = copy.State;
	N = State->GetParticleNumber();
	OccupancyList = new int[N];
	for(int i = 0 ; i<N ; i++)
		OccupancyList[i] = -1;
	CellList = NULL;
	DualList = NULL;
	return *this;
}

template <int Dim>
void CGrid<Dim>::SetState(CStaticState<Dim> *s)
{
	State = s;
}


//Functions to construct the grid
template <int Dim>
void CGrid<Dim>::Allocate()
{
	CellSize = ComputeCellSize();

	//Compute the number of grid elements in each direction. Additionally,
	//figure out the total number of cells as the product of all of these.
	int product = 1;
	for(int i = 0 ; i < Dim ; i++){
		N_Cells[i] = (int)ceil(1.0/CellSize(i));
		product*=N_Cells[i];
	}

	TotalCells = product;

	if(!(CellList==NULL))
		delete[] CellList;

	if(!(DualList==NULL))
		delete[] DualList;

	CellList = new int[product];
	for(int i = 0 ; i< product ;i++)
		CellList[i] = -1;
	
	DualList = new list<int>[product];
}

template <int Dim>
void CGrid<Dim>::Construct()
{
	//reset the lists
	for(int i = 0 ; i<N ; i++)
		OccupancyList[i] = -1;

	for(int i = 0 ; i< TotalCells ;i++)
	{
		CellList[i] = -1;
		DualList[i].erase(DualList[i].begin(),DualList[i].end());
	}
	

	//compute the transformed positions of the particles.
	Eigen::VectorXd Positions;
	State->GetPositionsVirtual(Positions);
	
	//First put the particles in the cells.
	for(int i = 0 ; i<N ; i++)
	{
		//Calculate the index of the cell that we should be in
		int Cell_Index = 0;
		int Prod = 1;
		for(int j = 0 ; j<Dim ; j++)
		{
			Cell_Index += (int)(Prod*floor(Positions(Dim*i+j)/CellSize(j)));
			Prod*=N_Cells[j];
		}
		
		//Add the particle to that cell
		if(CellList[Cell_Index]==-1)
		{
			CellList[Cell_Index] = i; 
		}else{
			int curr = CellList[Cell_Index];
			while(OccupancyList[curr]!=-1)
				curr = OccupancyList[curr];
			OccupancyList[curr] = i;
		}
	}

	//Now construct the dual list
	//go through each pair of cells and find their coordinates
	dvec Displacement;
	for(int i = 0 ; i < TotalCells ; i++)
	{
		dvec CoordinateI;
		CellToCoordinates(i,CoordinateI);
		for(int j = i ; j < TotalCells ; j++)
		{
			dvec CoordinateJ;
			CellToCoordinates(j,CoordinateJ);
			
			State->GetBox()->MinimumDisplacement(CoordinateI,CoordinateJ,Displacement);
			
			//go through the dimensions and if any dimension is closer than CellSize+\epsilon apart
			//add the cells to each others dual
			bool add = true;
			for(int k = 0 ; k < Dim ; k++)
				if(abs(Displacement(k))>CellSize(k)+CellSize(k)/100.0)	
					add = false;
	
			if(add)
			{
				//cout << "Displacement = " << max_displacement << " : CellSize = " << CellSize(0) << endl;
				DualList[i].push_back(j);
				if(j!=i)
					DualList[j].push_back(i);
			}
			
		}

	}

}

//Functions to access the grid
template<int Dim>
void CGrid<Dim>::CellToCoordinates(int i,dvec &Coordinate)
{
	int Prod = 1;
	int LProd = 1;
	dbl Sub = 0;
	for(int k = 0 ; k < Dim ; k++){
		Prod*=N_Cells[k];
		Coordinate(k) = (i%Prod) - Sub;
		Sub+=Coordinate(k);
		Coordinate(k)/=LProd;
		LProd*=N_Cells[k];
	}

	for(int i = 0; i< Dim ; i++)
		Coordinate(i)*=CellSize(i);
}
	
template <int Dim>
int CGrid<Dim>::CoordinatesToCell(const dvec &coordinates)
{
	int Cell_Index = 0;
	int Prod = 1;
	for(int j = 0 ; j<Dim ; j++)
	{
		Cell_Index += (int)(Prod*floor(coordinates(j)/CellSize(j)));
		Prod*=(int)N_Cells[j];
	}
	
	return Cell_Index;
}

template <int Dim>
typename CGrid<Dim>::iterator CGrid<Dim>::GetParticleIterator(int i)
{
	dvec pos;
	State->GetParticlePositionVirtual(pos,i);
	int j = CoordinatesToCell(pos);
	typename CGrid<Dim>::iterator ret(DualList[j],OccupancyList,CellList);
	return ret;
}

template <int Dim>
Eigen::Matrix<dbl,Dim,1> CGrid<Dim>::ComputeCellSize()
{
	dbl MaximumRadius = 0.0;
	Eigen::VectorXd Rads;
	State->GetRadii(Rads);

	for(int i = 0; i < State->GetParticleNumber() ; i++)
		if(Rads(i)>MaximumRadius)
			MaximumRadius = Rads(i);

	dbl Scale = 2*MaximumRadius;
	dbl N = floor(1.0/Scale);
	Scale = 1.0/N;

	dvec ret = dvec::Constant(Scale);
	
	return ret;
}

}

#endif
