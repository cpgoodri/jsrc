#include <iostream>
#include "State/StaticState.h"
#include "Potentials/RegisterPotentials.h"
#include "Boundaries/RegisterBoxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
using namespace std;

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
	pos->get(pv.data(),1,N,2);
	
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

int main()
{
	LiuJamming::RegisterPotentials();
	LiuJamming::RegisterBoxes<2>();

	NcError error(NcError::silent_nonfatal);

	LiuJamming::CStaticState<2> System(32);

	NcFile file("/data1/jcode/Carl_2012_states/4set-2d_harmonic_poly_uniform/state_bin/statedb_N00032_Lp-1.0000.nc");

	LoadCarlSystem(file,0,System);

	cout << "System created.\n";
	cout << "Particle positions:\n";
	Eigen::VectorXd Positions;
	System.GetPositions(Positions);
	cout << Positions << endl;
	
	LiuJamming::CStaticComputer<2> Computer(System);
	cout << "Energy = " << Computer.ComputeEnergy() << endl;
	
	
	return 0;
}
