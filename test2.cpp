#include <iostream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
//#include "Potentials/RegisterPotentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
#include <netcdf>

using namespace LiuJamming;
#define DIM 2
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

void test1(int N, dbl phi, int seed)
{
/*
	dbl a1=1;
	string s = ConvertDblToHexString(a1);
	cout << "s  = " << s << endl;

	dbl a2 = ConvertHexStringToDbl(s);
	cout << "a1 = " << a1 << endl;
	cout << "a2 = " << a2 << endl;
	cout << "a2-a1 = " << a2-a1 << endl;
*/

	/*
	cout << "Pot test:" << endl;
	CHarmonicPotential pot(1.2134);
	string str = pot.DataToString();
	cout << "str = " << str << endl;

	cout << "\nBox test:" << endl;
	CPeriodicBox<DIM> box;
	box.SetSymmetry(box.RECTANGULAR);
	str = box.DataToString();
	cout << "str = " << str << endl;
*/

	//Create the system
	CStaticState<DIM> s(N);
	s.RandomizePositions(seed);
	s.SetRadiiPolyUniform();
	s.SetPackingFraction(phi);

	CStaticState<DIM> s2 = s;


	CStaticComputer<DIM> c(s);
	CSimpleMinimizer<DIM> miner(c, CSimpleMinimizer<DIM>::FIRE);
	
//	CStaticComputer<DIM> c2(s2);
//	CSimpleMinimizer<DIM> miner2(c2, CSimpleMinimizer<DIM>::FIRE);

try{
	string fn = "test.nc";
	NcFile outFile(fn, NcFile::replace);
	NcDim recDim = outFile.addDim("rec");
	NcVar myVar = outFile.addVar("testVar", ncDouble, recDim);
}
catch (NcException& e)
{
	cout << "unknown error"<<endl;
	e.what();
}


	/*
	//Create the system
	CStaticState<DIM> System(N);

	//Set the initial positions... either according to some lattice or random.
	System.RandomizePositions(seed);
//	System.SetSquareLattice();

	//Set the radii distribution.
	System.SetRadiiPolyUniform();

	//Set the packing fraction.
	System.SetPackingFraction(phi); //This only changes the box size!

	//////////////////////////////////////////////
	//The system is now created and initialized.//
	//////////////////////////////////////////////

	//Print the volume
	printf("Volume = %e\n", System.GetVolume());

	//Minimize the energy:
	//	  First, create a Computer object
	CStaticComputer<DIM> Computer(System);
	//	  Then Create a Minimizer object
	CSimpleMinimizer<DIM> miner(Computer);
	//	  Then call the FIRE minimization routine.
	miner.minimizeFIRE();
	//miner.minimizeFIRE(1e-12,-1,-1.,500);
	//Alternatively, the above two lines can be combined into one:
	//CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);


	*/
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



	test1(N, phi, r);


	



}








