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

using namespace LiuJamming;
#define DIM 2
using namespace std;

typedef Eigen::Matrix<dbl,DIM,DIM> dmat;

/*
void test1(int N, dbl phi, int seed)
{




	//db.OpenReplace();

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
*/

void test2(int N, dbl phi, int seed)
{
	//Create the system
	CStaticState<DIM> s(N);

	s.RandomizePositions(seed);
	s.SetRadiiPolyUniform();
	s.SetPackingFraction(phi);

	CStaticComputer<DIM> c(s);
	CSimpleMinimizer<DIM> miner(c, CSimpleMinimizer<DIM>::FIRE);

	//Construct a bonds list and compute the bonds.
	CBondList<DIM> bonds;
	c.ComputeBondList(bonds);
	bonds.RemoveRattlers(DIM+1,true);

	MatrixInterface<dbl> MI;
	dmat strain = 0.5*dmat::Identity();
	bonds.ComputeGeneralizedHessian(MI.A, strain);

	MI.VDiagonalize();


}


int main(int argc, char* argv[])
{
	int N = 256;			//n
	dbl Lp = -1.0;			//p
	dbl phi = 0.9;			//f
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








