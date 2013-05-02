#include <iostream>
#include "State/StaticState.h"
#include "Potentials/RegisterPotentials.h"
#include "Boundaries/RegisterBoxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
using namespace std;
using namespace LiuJamming;
#define DIM 2

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
	LiuJamming::RegisterPotentials();
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
	printf("number of bonds = %i\n", bonds.GetNBonds());
	cout << "Energy = " << bonds.ComputeEnergy() << endl;

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
	
	return 0;
}

/*
void test1()
{
	RegisterPotentials();
	RegisterBoxes<DIM>();

	int N=32;
	long seed = 123;
	CStaticState<DIM> System(N);
	System.RandomizePositions();

	System.SetPackingFraction(0.9);
	System.PrintParticles();
	printf("Volume = %e\n", System.GetVolume());

	CStaticComputer<DIM> Computer(System);
	
	CBondList<DIM> bonds;
	Computer.ComputeBondList(bonds);
	printf("number of bonds = %i\n", (int)bonds.list.size());
	bonds.PrintBonds();
	
//	CSimpleMinimizer<DIM> miner(Computer);
//	miner.minimizeFIRE(1e-12,0.1);
}
*/










int main(int argc, char* argv[])
{
	int N = 256;			//n
	dbl Lp = -1.0;			//p
	int r = 1;				//r

	int c;
	while((c=getopt(argc, argv, "n:p:r:")) != -1)
		switch(c)
		{
			case 'n':	N = atoi(optarg); break;
			case 'p':	Lp = atof(optarg); break;
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


	//test1();
	SamTest(N, Lp, r);


	



}








