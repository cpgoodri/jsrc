//Constructor and copy operators
template <int Dim>
CGrid<Dim>::CGrid(int _N, const VectorXd &_ParticlePositions, CBox *_Box) : N(_N), ParticlePositions(_ParticlePositions), Box(_Box)
{
	OccupancyList = new int[N];
	for(int i = 0 ; i<N ; i++)
		OccupancyList[i] = -1;
	CellList = NULL;
	DualList = NULL;
}

template <int Dim>
CGrid<Dim>::CGrid(const CGrid &copy)
{

}

template <int Dim>
const CGrid<Dim>::CGrid &operator=(const CGrid &copy)
{

}
	
//Functions to set the various properties of the grid
template <int Dim>
void CGrid<Dim>::SetCellSize(double size)
{
	CellSize = size;
}

//Functions to construct the grid
template <int Dim>
void CGrid<Dim>::AllocateGrid()
{
	//Compute the number of grid elements in each direction. Additionally,
	//figure out the total number of cells as the product of all of these.
	Eigen::Matrix<double,Dim,1> corner = Eigen::Matrix<double,Dim,1>::Constant(1.0);
	Box->Transform(corner);
	double product = 1;
	for(int i = 0 ; i < Dim ; i++){
		N_Cells[i] = ceil(corner(i)/CellSize);
		prod*=N_Cells[i];
	}

	TotalCells = prod;

	CellList = new int[prod];
	for(int i = 0 ; i< prod ;i++)
		CellList[i] = -1;
	DualList = new list<int>[prod];
}

template <int Dim>
void CGrid<Dim>::ConstructGrid()
{
	//reset the lists
	for(int i = 0 ; i<N ; i++)
		OccupancyList[i] = -1;

	for(int i = 0 ; i< TotalCells ;i++)
	{
		CellList[i] = -1;
		DualList.erase(DualList.begin(),DualList.end());
	}
	

	//compute the transformed positions of the particles.
	VectorXd T_Positions = ParticlePositions;
	Box->Transform(T_Positions);
	
	//First put the particles in the cells.
	for(int i = 0 ; i<N ; i++)
	{
		//Calculate the index of the cell that we should be in
		int Cell_Index = 0;
		int Prod = 1;
		for(int j = 0 ; j<Dim ; j++)
		{
			Cell_Index += Prod*floor(T_Positions(Dim*i+j)/CellSize);
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
	Eigen::Matrix<double,Dim,1> Displacement;
	for(int i = 0 ; i < TotalCells ; i++)
	{
		Eigen::Matrix<double,Dim,1> CoordinateI;
		CellToCoordinates(i,CoordinateI);

		for(int j = i ; j < TotalCells ; j++)
		{
			Eigen::Matrix<double,Dim,1> CoordinateJ;
			CellToCoordinates(i,CoordinateJ);
			
			Box->MinimumDisplacement(CoordinateI,CoordinateJ,Displacement);
			
			//go through the dimensions and if any dimension is closer than CellSize+\epsilon apart
			//add the cells to each others dual
			for(int k = 0 ; k < Dim ; k++)
				if(abs(Displacement)<=CellSize+CellSize/100.0)
				{
					DualList[i].push_back(j);
					if(j!=i)
						DualList[j].push_back(i);
				}
			
		}
	}
	
	
}

//Functions to access the grid
template<int Dim>
const Matrix<double,Dim,1> &CGrid::CellToCoordinates(int i,Eigen::Matrix<double,Dim,1> &Coordinate)
{
	int Prod = 1;
	int LProd = 1;
	int Sub = 0;
	for(int k = 0 ; k < Dim ; k++){
		Prod*=N_Cells[j];
		CoordinateI(k) = (i%Prod) - Sub;
		Sub+=CoordinateI(k);
		CoordinateI(k)/=LProd;
		LProd*=N_Cells[j];
	}
	
	Coordinate*=CellSize;
	
	return Coordinate;
}
	
template <int Dim>
int CGrid::CoordinatesToCell(const Matrix<double,Dim,1> &coordinates)
{
	int Cell_Index = 0;
	int Prod = 1;
	for(int j = 0 ; j<Dim ; j++)
	{
		Cell_Index += Prod*floor(coordinates(j)/CellSize);
		Prod*=N_Cells[j];
	}
	
	return Cell_Index;
}
	
	
	