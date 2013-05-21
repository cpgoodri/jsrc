#include <iostream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
//#include "Potentials/RegisterPotentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
using namespace std;
using namespace LiuJamming;
#define DIM 2



/*
double GetScale(int N,double phi, Eigen::VectorXd &rv)
{
	assert(DIM==2);
	dbl phi0 = M_PI*rv.dot(rv);
	dbl r0 = sqrt(phi0/phi);
	return 1./r0;
//	return sqrt(phi/(M_PI*N*0.5*(1+1.4*1.4)));
}

void LoadCarlSystem(NcFile &file,int Record,LiuJamming::CStaticState<2> &sys)
{
	NcVar *pos = file.get_var("pos");
	NcVar *radii = file.get_var("relRad");
	NcVar *pf = file.get_var("pf");
	NcVar *energy = file.get_var("energy");

	int N = file.get_dim("NP")->size();

	Eigen::VectorXd pv = Eigen::ArrayXd::Zero(DIM*N);
	Eigen::VectorXd rv = Eigen::ArrayXd::Zero(N);
	double p = 0.0;
	double en = 0.0;

	pos->set_cur(Record);
	pos->get(pv.data(),1,N*DIM);
	
	radii->get(rv.data(),N);

	pf->set_cur(Record);
	pf->get(&p,1);

	energy->set_cur(Record);
	energy->get(&en,1);

	cout << "System at packing fraction " << p << " with energy " << en << endl;

	sys.SetPositions(pv);
	dbl scale = GetScale(N,p,rv);
	rv*=scale;
	sys.SetRadii(rv);

//	sys.SetPackingFraction(p);

}

int SamTest(int N=32, dbl Lp=-1.0, int rec=0)
{
	//LiuJamming::RegisterPotentials();
	LiuJamming::RegisterBoxes<DIM>();

	NcError error(NcError::silent_nonfatal);

	LiuJamming::CStaticState<DIM> System(N);

	char filename[256];
	sprintf(filename,"/data1/jcode/Carl_2012_states/4set-2d_harmonic_poly_uniform/state_bin/statedb_N%05i_Lp%6.4f.nc", N, Lp);
	NcFile file(filename);		
	//NcFile file("/data1/jcode/Carl_2012_states/4set-2d_harmonic_poly_uniform/state_bin/statedb_N00032_Lp-1.0000.nc");

	LoadCarlSystem(file,rec,System);

	cout << "System created.\n";
//	cout << "Particle positions:\n";
//	Eigen::VectorXd Positions;
//	System.GetPositions(Positions);
//	cout << Positions << endl;

//	System.PrintParticles();
	printf("Volume = %e\n", System.GetVolume());
	
	LiuJamming::CStaticComputer<DIM> Computer(System);
//	cout << "Energy = " << Computer.ComputeEnergy() << endl;

	CBondList<2> bonds;
	Computer.ComputeBondList_NoGrid(bonds);
	printf("number of bonds  = %i\n", bonds.GetNBonds());
	cout << "Energy  = " << bonds.ComputeEnergy() << endl;

	fflush(stdout);
	CBondList<2> bonds2;
	Computer.ComputeBondList(bonds2);
	printf("number of bonds2 = %i\n", bonds2.GetNBonds());
	cout << "Energy2 = " << bonds2.ComputeEnergy() << endl;

	if(false)
	{
		bonds.RemoveRattlers(DIM+1,true);
		printf("number of bonds = %i\n", bonds.GetNBonds());
		cout << "Energy = " << bonds.ComputeEnergy() << endl;

		Eigen::SparseMatrix<dbl> Smat;
		bonds.ComputeHessian(Smat);


	//	Eigen::VectorXd grad;
	//	dbl energy = bonds.ComputeGradient(grad);
	//	cout << "grad:\n" << grad << endl;

		MatrixInterface<dbl> mat;
		bonds.ComputeHessian(mat.A);
		mat.Diagonalize();
		//mat.HumanReport();
		for(int i=0; i<mat.num_converged; ++i)
			printf(" eigenvalue %5i = % e\n", i, mat.Eigenvalues[i]);
	}

	return 0;
}
*/

void test1(int N, dbl phi, int seed)
{
	//You still need to register boxes... This will change.
//	RegisterBoxes<DIM>();

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

	//Print the particles
	System.PrintParticles();

	//Print the volume
	printf("Volume = %e\n", System.GetVolume());

	//Minimize the energy:
	//	  First, create a Computer object
	CStaticComputer<DIM> Computer(System);
	//	  Then Create a Minimizer object
	CSimpleMinimizer<DIM> miner(Computer);
	//	  Then call the FIRE minimization routine.
	miner.minimizeFIRE();
	//Alternatively, the above two lines can be combined into one:
	//CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);




	//////////////////////////////////////////////////////////////////////////////
	//Now do calculations. Much of this will be moved into the Computer class...//
	//////////////////////////////////////////////////////////////////////////////

	//Construct a bonds list and compute the bonds.
	CBondList<DIM> bonds;
	Computer.ComputeBondList(bonds);
	
	//Remove rattlers
	bonds.RemoveRattlers(DIM+1,true);

	//Simple calcualtions can be made directly to the bonds list:
	printf("number of bonds = %i\n", bonds.GetNBonds());
	printf("pressure = %e\n", bonds.ComputePressure());

	//If uncommented, the following line would remove the stress from the bonds list.
//	bonds.MakeUnstressed();

	//Calculate and print the elastic constants
	cCIJKL<DIM> cijkl;
	bonds.CalculateCijkl(cijkl);
	printf("bulk  modulus = % e\n", cijkl.CalculateBulkModulus());
	printf("shear modulus = % e\n", cijkl.CalculateShearModulus());
	printf("weakest deformation = % e\n", cijkl.calc_weakest_deformation());


	//To find eigenmodes, first construct a MatrixInterface
	MatrixInterface<dbl> D1;
	
	//Set the hessian from the bonds list
	bonds.ComputeHessian(D1.A);

	//Diagonalize the hessian
	D1.Diagonalize();			//Get all Eigenvalues/vectors (except the largest one)
//	D1.VDiagonalize();			//Same as above, but print Eigenvalue Report at end
//	D1.Diagonalize(100);		//Get only the lowest 100
//	D1.VDiagonalize(100);		//Same as above, but print Eigenvalue Report at end


	//Report the first 8 (default=30) eigenvalues (plus some header info) to stdout.
	D1.Report(8);
}


void test2(int N, dbl phi, int seed)
{
	//You still need to register boxes... This will change.
//	RegisterBoxes<DIM>();

	//Create the system
	CStaticState<DIM> System(N);
	System.RandomizePositions(seed);
	System.SetRadiiPolyUniform();
	System.SetPackingFraction(phi);

	//Minimize the energy:
	CStaticComputer<DIM> Computer(System);
	CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);

	//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	Computer.StdPrepareSystem();

	//Calculate info
	printf("\nCalculate data for stressed system-->\n");
	Computer.CalculateStdData(); //Have option to not write down hessian.
	Computer.Data.Print();
	Computer.Data.H.VDiagonalize();

	//Calculate unstressed info
	printf("\nCalculate data for unstressed system-->\n");
	CStdData<DIM> unstressedData;
	Computer.CalculateStdData_Unstressed(unstressedData);
	unstressedData.Print();
	unstressedData.H.VDiagonalize();
}



void test3(int N, dbl phi, int seed)
{
	//You still need to register boxes... This will change.
//	RegisterBoxes<DIM>();

	//Create the system
	CStaticState<DIM> System(N);
	System.RandomizePositions(seed);
	System.SetRadiiPolyUniform();
	System.SetPackingFraction(phi);

	//Minimize the energy:
	CStaticComputer<DIM> Computer(System);
	CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);

	//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	Computer.StdPrepareSystem();

	//Diagonalize H(k)
	MatrixInterface<cdbl> Dk;
	Eigen::SparseMatrix<cdbl> Tk;
	typedef Eigen::Matrix<dbl, DIM, 1> dvec;
	dvec k = dvec::Zero();
	k[0] = 1e-2;
	Computer.ComputeHessianBZ(Dk.A, Tk, k);
	Dk.VDiagonalize();
	cout << Dk.A << endl;
}



void test4(int seed)
{
	int N = 512;
	dbl phi = 0.9; 

	//You still need to register boxes... This will change.
//	RegisterBoxes<DIM>();

	//Create the system
	CStaticState<DIM> System(N);
	System.RandomizePositions(seed);
	System.SetRadiiBi();
	System.SetPackingFraction(phi);

	if(true)
	{
		//Minimize the energy:
		CStaticComputer<DIM> Computer(System);
		CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);

		//Construct a bonds list and compute the bonds.
		CBondList<DIM> bonds;
		Computer.ComputeBondList(bonds);
		
		//Remove rattlers
		bonds.RemoveRattlers(DIM+1,true);

		//Simple calcualtions can be made directly to the bonds list:
		printf("number of bonds = %i\n", bonds.GetNBonds());
		printf("pressure = %e\n", bonds.ComputePressure());

		char filename[256];
		sprintf(filename, "pos_test.txt");
		FILE *outfile1 = fopen(filename, "w");
		
		Eigen::VectorXd pos, rad;
		System.GetPositions(pos);
		System.GetRadii(rad);
		for(int i=0; i<System.GetParticleNumber(); ++i)
			fprintf(outfile1, "% e % e % e\n", pos[DIM*i], pos[DIM*i+1], rad[i]);
		fflush(outfile1);
		fclose(outfile1);

		sprintf(filename, "bonds_test.txt");
		bonds.PrintBonds_txt(filename);
	}
	if(true)
	{
		System.SetPackingFraction(0.8);

		//Minimize the energy:
		CStaticComputer<DIM> Computer(System);
		CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);
		
		//Construct a bonds list and compute the bonds.
		CBondList<DIM> bonds;
		Computer.ComputeBondList(bonds);
		
		//Remove rattlers
//		bonds.RemoveRattlers(DIM+1,true);

		char filename[256];
		sprintf(filename, "pos_test2.txt");
		FILE *outfile2 = fopen(filename, "w");
		
		Eigen::VectorXd pos, rad;
		System.GetPositions(pos);
		System.GetRadii(rad);
		for(int i=0; i<System.GetParticleNumber(); ++i)
			fprintf(outfile2, "% e % e % e\n", pos[DIM*i], pos[DIM*i+1], rad[i]);
		fflush(outfile2);
		fclose(outfile2);



	}

}






int main(int argc, char* argv[])
{
	int N = 256;			//n
	dbl Lp = -1.0;			//p
	dbl phi = 0.9;			//f
	int r = 1;				//r

	int c;
	while((c=getopt(argc, argv, "n:p:f:r:")) != -1)
		switch(c)
		{
			case 'n':	N = atoi(optarg); break;
			case 'p':	Lp = atof(optarg); break;
			case 'f':	phi = atof(optarg); break;
			case 'r':	r = atoi(optarg); break;
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


//	test1(N, phi, r);
	test2(N, phi, r);
//	test3(N, phi, r);
//	test4(r);
	//SamTest(N, Lp, r);


	



}








