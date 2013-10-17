#include <iostream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
//#include "Potentials/RegisterPotentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
#include "Database/Database.h"

#include "Resources/Clusterizer.h"

using namespace LiuJamming;
#define DIM 2
using namespace std;

typedef Eigen::Matrix<dbl,DIM,DIM> dmat;


void test1(int N, dbl phi, int seed)
{
	//Create the system
	CStaticState<DIM> s(N);
if(false)
{
	s.RandomizePositions(seed);
	s.SetRadiiPolyUniform();
	s.SetPackingFraction(phi);

	CStaticDatabase<DIM> db(N,"temp2.nc",NcFile::Replace);
	db.WriteState(s);

	CStaticComputer<DIM> c(s);
	CSimpleMinimizer<DIM> miner(c, CSimpleMinimizer<DIM>::FIRE);
	
	db.WriteState(s, 4);
}else{
	CStaticDatabase<DIM> db(N,"temp2.nc",NcFile::ReadOnly);
	db.ReadState(s,0);
	
	CStaticComputer<DIM> c(s);
	CSimpleMinimizer<DIM> miner(c, CSimpleMinimizer<DIM>::FIRE);
}
}







void test2(int N, dbl phi, int seed)
{
	//Create the system
	CStaticState<DIM> s(N);

	s.RandomizePositions(seed);
	s.SetRadiiPolyUniform();
	s.SetPackingFraction(phi);

	CStaticComputer<DIM> c(s);
	CSimpleMinimizer<DIM> miner(c, CSimpleMinimizer<DIM>::FIRE);
	
	//Prepare system for standard calculations
	//CStaticComputer<DIM>::StdPrepareSystem() sets the internal bond list and removes rattlers
	//it returns 0 if everything is done correctly.
	if(c.StdPrepareSystem())
	{
		printf("The bonds list is empty...\n");
		return;
	}
	c.CalculateStdData(); //Have option to not write down hessian.
	c.Data.Print();


	vector< vector<int> > nbrs;
	c.Bonds.CalculateNeighbors(nbrs);

	//randomly set some nodes to "include"
	std::vector<bool> include(c.Bonds.GetN(), false);
	for(int i=0; i<c.Bonds.GetN()/2; ++i)
		include[i] = true;

	//Create a clusterizer and run the decomposition
	CClusterizer<DIM> cizer(c.Bonds.GetN(), nbrs, include);
	printf("Number of clusters = %i\n", cizer.GetNumClusters());

	//Get all the cluseters
	vector< vector<int> > clusters;
	cizer.GetClusters(clusters);

	//Print the clusters
	for(int c=0; c<(int)clusters.size(); ++c)
	{
		printf(" cluster %5i: ", c);
		for(int i=0; i<(int)clusters[c].size(); ++i)
			printf("%i ", clusters[c][i]);
		printf("\n");
	}


	//Construct a bonds list and compute the bonds.
//	CBondList<DIM> bonds;
//	c.ComputeBondList(bonds);
//	bonds.RemoveRattlers(DIM+1,true);

	/*
	MatrixInterface<dbl> MI;
	dmat strain = dmat::Zero();
	strain = 0.5*dmat::Identity();
	c.Bonds.ComputeGeneralizedHessian(MI.A, strain);

	MI.VDiagonalize_Eigen();
	int Nvar = DIM*c.Data.NPp+1;
	int nc = MI.num_converged;

	for(int m=0; m<50; ++m)
	{
		printf("%5i   % e   % e\n", m, MI.Eigenvalues[m], MI.Eigenvalues[m]/POW2(MI.Eigenvectors[Nvar*m+(Nvar-1)]));
	}
	*/

	/*
	dbl G = 0.;
	assert(nc == Nvar);
	for(int m=0; m<nc; ++m)
	{
		G += POW2(MI.Eigenvectors[Nvar*m+(Nvar-1)])/MI.Eigenvalues[m];
	}
	G = 1./G;
	G /= c.GetVolume();

	dbl cxyxy = c.Data.cijkl.cxyxy;
	printf("G     = % e\n", G);
	printf("cxyxy = % e\n", cxyxy);
	printf("diff  = % e\n", G-cxyxy);

	*/



















}


int main(int argc, char* argv[])
{
	int N = 256;			//n
	dbl Lp = -1.0;			//p
	dbl phi = 0.845;			//f
	int r = 1;				//r
	int NFixedParticles = 0;//x

	int c;
	while((c=getopt(argc, argv, "n:p:f:r:x:")) != -1)
		switch(c)
		{
			case 'n':	N = atoi(optarg); break;
			case 'p':	Lp = atof(optarg); break;
			case 'f':	phi = atof(optarg); break;
			case 'r':	r = atoi(optarg); break;
			case 'x':	NFixedParticles = atoi(optarg); break;
			case '?':	
				if(optopt == 'c') 
					std::cerr << "Option -" << optopt << "requires an argument.\n";
				else if(isprint(optopt)) 
					std::cerr << "Unknown opton '-" << optopt << "'.\n";
				else
					std::cerr << "Unknown option character '\\x" << optopt << "'.\n";
				return 1;
			default:
				abort();
		}
	for(int index = optind; index < argc; index++)
		std::cerr << "Non-option argument " << argv[index] << "\n";



	test2(N, phi, r);


	



}








