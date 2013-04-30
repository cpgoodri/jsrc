#include <iostream>
#include "State/StaticState.h"
#include "Potentials/RegisterPotentials.h"
#include "Boundaries/RegisterBoxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Minimization/minimizer.h"
using namespace std;
using namespace LiuJamming;
#define DIM 2

double GetScale(int N,double phi)
{
	return sqrt(phi/(M_PI*N*0.5*(1+1.4*1.4)));
}

void LoadCarlSystem(NcFile &file,int Record,LiuJamming::CStaticState<2> &sys)
{
	NcVar *pos = file.get_var("pos");
	NcVar *radii = file.get_var("relRad");
	NcVar *pf = file.get_var("pf");
	NcVar *energy = file.get_var("energy");

	int N = file.get_dim("NP")->size();

	Eigen::VectorXd pv = Eigen::ArrayXd::Zero(2*N);
	Eigen::VectorXd rv = Eigen::ArrayXd::Zero(N);
	double p = 0.0;
	double en = 0.0;

	pos->set_cur(Record);
	pos->get(pv.data(),1,N*2);
	
	radii->get(rv.data(),N);

	pf->set_cur(Record);
	pf->get(&p,1);

	energy->set_cur(Record);
	energy->get(&en,1);

	cout << "System at packing fraction " << p << " with energy " << en << endl;

	sys.SetPositions(pv);
	rv*=GetScale(N,p);
	sys.SetRadii(rv);
	

}

int SamTest()
{
	LiuJamming::RegisterPotentials();
	LiuJamming::RegisterBoxes<2>();

	NcError error(NcError::silent_nonfatal);

	LiuJamming::CStaticState<2> System(32);

	NcFile file("/data1/jcode/Carl_2012_states/4set-2d_harmonic_poly_uniform/state_bin/statedb_N00032_Lp-1.0000.nc");

	LoadCarlSystem(file,0,System);

	cout << "System created.\n";
//	cout << "Particle positions:\n";
//	Eigen::VectorXd Positions;
//	System.GetPositions(Positions);
//	cout << Positions << endl;

	System.PrintParticles();
	printf("Volume = %e\n", System.GetVolume());
	
	LiuJamming::CStaticComputer<2> Computer(System);
	cout << "Energy = " << Computer.ComputeEnergy() << endl;

	CBondList<2> bonds;
	Computer.ComputeBondList(bonds);
	printf("number of bonds = %i\n", (int)bonds.list.size());
	
	
	bonds.PrintBonds();
	
	
	return 0;
}

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











int main()
{
	//test1();
	SamTest();


	



}








