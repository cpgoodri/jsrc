#include <iostream>
#include "State/StaticState.h"
#include "Potentials/RegisterPotentials.h"
#include "Boundaries/RegisterBoxes.h"
#include "Resources/Exception.h"
using namespace std;

int main()
{
	LiuJamming::RegisterPotentials();
	LiuJamming::RegisterBoxes<2>();

	NcError error(NcError::silent_nonfatal);

	LiuJamming::CStaticState<2> System(10);
	System.RandomizePositions();
	cout << "System created.\n";
	cout << "Particle positions:\n";
	Eigen::VectorXd Positions;
	System.GetPositions(Positions);
	cout << Positions << endl;
	cout << "In a sheared container the positions will be:\n";
	Eigen::Matrix<double,2,2> trans;
	trans(0,0) = 1;
	trans(1,1) = 1;
	trans(1,0) = 0;
	trans(0,1) = 1;
	System.GetBox()->SetTransformation(trans);
	cout << "New transformation matrix is:\n";
	Eigen::Matrix<double,2,2> trans2;
	System.GetBox()->GetTransformation(trans2);
	cout << trans2 << endl;
	cout << "New positions are:\n";
	System.GetPositions(Positions);
	cout << Positions << endl;
	
	cout << "Saving file\n";
	NcFile file("test.nc",NcFile::Replace);
	
	try{
		System.Write(file,0);
	}catch(LiuJamming::CException e){
		cout << e << endl;
	}


	cout << "Reading file.\n";
	NcFile file2("test.nc",NcFile::ReadOnly);
	
	LiuJamming::CStaticState<2> System2(10);
	try{
		System2.Read(file,0);
	}catch(LiuJamming::CException e){
		cout << e << endl;
	}
	
	cout << "Read system positions are:\n";
	System2.GetPositions(Positions);
	cout << Positions << endl;
	
	return 0;
}
