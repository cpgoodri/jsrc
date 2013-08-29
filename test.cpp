#include <iostream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
//#include "Potentials/RegisterPotentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
#include "Resources/OrderParams.h"


#include <Eigen/Eigenvalues>
using namespace std;
using namespace LiuJamming;
#define DIM 3



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

void test_orderparam(int N, dbl phi, int seed)
{
	//Create the system
	CStaticState<DIM> System(N);
	System.RandomizePositions(seed);
//	System.SetFCCLattice();
	System.SetRadiiMono();
	System.SetPackingFraction(phi);

	//Minimize the energy:
	System.AssumeRectangularSymmetry(); //This gives a slight performance boost.
	CStaticComputer<DIM> Computer(System);
	CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);
	
	//Prepare system for standard calculations
	//CStaticComputer<DIM>::StdPrepareSystem() sets the internal bond list and removes rattlers
	//it returns 0 if everything is done correctly.
	if(Computer.StdPrepareSystem())
	{
		printf("The bonds list is empty...\n");
		return;
	}

	//Calculate info
	printf("\nCalculate data for stressed system-->\n");
	Computer.CalculateStdData(); //Have option to not write down hessian.
	Computer.Data.Print();

	LocalOrder<DIM> lo(Computer.Bonds);
	lo.Calculate();

	cout << lo.fi << endl;


}

#if false

void test0(int N, dbl phi, int seed)
{
	//Create the system
	CStaticState<DIM> System(N);
	System.RandomizePositions(seed);
	System.SetRadiiPolyUniform();
	System.SetPackingFraction(phi);

	//Minimize the energy:
	System.AssumeRectangularSymmetry(); //This gives a slight performance boost.
	CStaticComputer<DIM> Computer(System);
	CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);
}

void test1(int N, dbl phi, int seed)
{
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
	//miner.minimizeFIRE(1e-12,-1,-1.,500);
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
//	D1.Diagonalize();			//Get all Eigenvalues/vectors (except the largest one)
//	D1.VDiagonalize();			//Same as above, but print Eigenvalue Report at end
//	D1.Diagonalize(100);		//Get only the lowest 100
//	D1.VDiagonalize(100);		//Same as above, but print Eigenvalue Report at end
	
	D1.VDiagonalize(100);		//Same as above, but print Eigenvalue Report at end


	//Report the first 8 (default=30) eigenvalues (plus some header info) to stdout.
	D1.Report(8);
}


void test2(int N, dbl phi, int seed)
{
	//Create the system
	CStaticState<DIM> System(N);
	System.RandomizePositions(seed);
	System.SetRadiiPolyUniform();
	System.SetPackingFraction(phi);

	//System.SetPotentialHertzian();
	//System.SetPotential(new CHertzianPotential(1.));
	System.SetPotential(new CSoftPotential(1.,2.5));

	//Minimize the energy:
	CStaticComputer<DIM> Computer(System);
	CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE);

//	System.SetBoxPeriodic(DIM);
//	System.SetPackingFraction(phi);
//	CStaticComputer<DIM> Computer2(System);

	//Prepare system for standard calculations
	//CStaticComputer<DIM>::StdPrepareSystem() sets the internal bond list and removes rattlers
	//it returns 0 if everything is done correctly.
	if(Computer.StdPrepareSystem())
	{
		printf("The bonds list is empty...\n");
		return;
	}

	//Calculate info
	printf("\nCalculate data for stressed system-->\n");
	Computer.CalculateStdData(); //Have option to not write down hessian.
	Computer.Data.Print();
	//Computer.Data.H.VDiagonalize();

	typedef Eigen::Matrix<dbl,DIM,DIM> dmat;
	dmat fab;
	Computer.Bonds.ComputeFabricTensor(fab);
	cout << "Fabric Tensor\n" << fab << endl;

	cout << "eigenvalues:\n" << fab.eigenvalues() << endl;
	cout << "Trace  = " << fab.trace() << endl;

	/*
	//Calculate unstressed info
	printf("\nCalculate data for unstressed system-->\n");
	CStdData<DIM> unstressedData;
	Computer.CalculateStdData_Unstressed(unstressedData);
	unstressedData.Print();
	unstressedData.H.VDiagonalize();
	*/
}


void test_fixed(int N, dbl phi, int seed, int Nfixed_particles = 0)
{
	//Create the system
	CStaticState<DIM> System(N);
	System.RandomizePositions(seed);
	System.SetRadiiPolyUniform();
	System.SetPackingFraction(phi);

	//You should actually fix random particles, not just the first ones.
	vector<bool> FixedDof(DIM*N,false); //1 -> fixed, 0 -> not fixed. size DIM*N
	for(int i=0; i<Nfixed_particles; ++i)
	{
		FixedDof[DIM*i  ] = true;
		FixedDof[DIM*i+1] = true;
	}

	//Minimize the energy:
	const dbl tol = 1e-12;
	CStaticComputer<DIM> Computer(System);
	Computer.SetFixedDof(FixedDof);
	CSimpleMinimizer<DIM> miner(Computer, CSimpleMinimizer<DIM>::FIRE, tol);

	//Prepare system for standard calculations 
	//	The routine to remove rattlers is currently only set up to deal with entirely fixed particles.
	//	We first convert FixedDof into FixedParticles and then pass this to the StdPrepareSystem routine
	vector<bool> FixedParticles(N,false);
	for(int i=0; i<N; ++i)
		if(FixedDof[DIM*i] || FixedDof[DIM*i+1])
			FixedParticles[i] = true;

	//CStaticComputer<DIM>::StdPrepareSystem() sets the internal bond list and removes rattlers
	//it returns 0 if everything is done correctly.
	if(Computer.StdPrepareSystem(FixedParticles))
	{
		printf("The bonds list is empty...\n");
		return;
	}

	//Calculate info
	printf("\nCalculate data for stressed system-->\n");
	Computer.CalculateStdData(false,false); //Don't compute elastic constants or write down the dynamical matrix.
	
	//NOTE: Here, the value of max_grad was calculated without knowing that particles are fixed. Redo this calculation manually:
	//	We first have to update FixedDof according to the rattler map (found in Computer.RattlerMap).
	assert(Computer.RattlerMap.full_size == N);
	vector<bool> FixedDofNew;
	for(int i=0; i<N; ++i)
		if(Computer.RattlerMap.inv(i) != -1)
		{
			FixedDofNew.push_back(FixedDof[DIM*i]);
			FixedDofNew.push_back(FixedDof[DIM*i+1]);
		}
	//	Now compute the gradient
	Eigen::VectorXd grad;
	Computer.Bonds.ComputeGradient(grad, FixedDofNew);
	Computer.Data.MaxGrad = max_abs_element(grad);


	//Print out data to the screen
	Computer.Data.Print();

	int Nc_mm=0, Nc_mf=0, Nc_ff=0;
	bool ifixed, jfixed;
	for(vector< CBond<DIM> >::const_iterator b=Computer.Bonds.begin(); b!=Computer.Bonds.end(); ++b)
	{
		assert(b->i >= 0 && b->i < Computer.RattlerMap.size());
		assert(b->j >= 0 && b->j < Computer.RattlerMap.size());
		ifixed = FixedParticles[Computer.RattlerMap[b->i]];
		jfixed = FixedParticles[Computer.RattlerMap[b->j]];

		if(ifixed && jfixed)
			++Nc_ff;
		else if(ifixed || jfixed)
			++Nc_mf;
		else
			++Nc_mm;
	}

	int min_passed = (Computer.Data.MaxGrad < tol)?1:0;

	//Save data. Right now this prints to the screen, but could easily have it print to a file.

	//DATA: N, phi, Nfixed, seed, energy, pressure, max_grad, min_passed, NPp, Nc, Nc_mm, Nc_mf, Nc_ff
	printf("%5i %f %5i %5i % e % e % e %5i %5i %5i %5i %5i %5i\n", N, phi, Nfixed_particles, seed, Computer.Data.Energy, Computer.Data.Pressure, Computer.Data.MaxGrad, min_passed, Computer.Data.NPp, Computer.Data.Nc, Nc_mm, Nc_mf, Nc_ff);
	fflush(stdout);
	//
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

	//Prepare system for standard calculations
	//CStaticComputer<DIM>::StdPrepareSystem() sets the internal bond list and removes rattlers
	//it returns 0 if everything is done correctly.
	if(Computer.StdPrepareSystem())
	{
		printf("The bonds list is empty...\n");
		return;
	}


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



#endif //false


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



//	test_fixed(N, phi, r, NFixedParticles);

	test_orderparam(N, phi, r);

//	test0(N, phi, r);
//	test1(N, phi, r);
//	test2(N, phi, r);
//	test3(N, phi, r);
//	test4(r);
	//SamTest(N, Lp, r);


	



}








