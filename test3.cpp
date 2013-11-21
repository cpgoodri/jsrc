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

#include "Resources/MatrixMarket.h"
#include "Resources/Clusterizer.h"

using namespace LiuJamming;
#define DIM 2
using namespace std;

typedef Eigen::Matrix<dbl,DIM,DIM> dmat;

void test2(int N, dbl phi, int seed)
{
	//Create the system
	CStaticState<DIM> s(N);
	s.RandomizePositions(seed);
	s.SetRadiiMono();
	s.SetPackingFraction(phi);
	//s.SetHexLattice();

	CStaticComputer<DIM> c(s);
	CSimpleMinimizer<DIM> miner(c, CSimpleMinimizer<DIM>::FIRE);
	if(c.StdPrepareSystem())
	{
		printf("The bonds list is empty...\n");
		return;
	}
	c.CalculateStdData(false);
	c.Data.Print();

	Eigen::Matrix<int,DIM,1> copies = Eigen::Matrix<int,DIM,1>::Constant(10);
	CStaticState<DIM> s2(s, copies);

	CStaticComputer<DIM> c2(s2);
	if(c2.StdPrepareSystem())
	{
		printf("The bonds list is empty...\n");
		return;
	}
	c2.CalculateStdData(false);
	c2.Data.Print();




	MatrixInterface<dbl> H0;
	c2.Bonds.ComputeHessian(H0.A, 1., 1., 1e-12);
	MatrixMarket::WriteSparseDbl("scratch/MMtest_H0.mtx", H0.A, true);

	MatrixInterface<dbl> Hg;
	dmat strain_tensor = dmat::Zero();
	strain_tensor(0,1) = strain_tensor(1,0) = 1.;
	c2.Bonds.ComputeGeneralizedHessian(Hg.A, strain_tensor);
	assert(Hg.A.rows()-1 == H0.A.rows());

	Eigen::VectorXd temp = Eigen::VectorXd::Zero(Hg.A.rows());
	temp[Hg.A.rows()-1] = 1.;

	Eigen::VectorXd temp2 = Hg.A*temp;
	Eigen::VectorXd F = Eigen::VectorXd::Zero(Hg.A.rows()-1);
	for(int i=0; i<Hg.A.rows()-1; ++i)
		F[i] = temp2[i];

	MatrixMarket::WriteDenseDbl("scratch/MMtest_F.mtx", F);

	H0.LUdecomp();
	Eigen::VectorXd X = Eigen::VectorXd::Zero(F.size());
	H0.solve_Mx_equals_b(X,F);

/*
	//Load tesla solution
	Eigen::VectorXd X_in;
	MatrixMarket::ReadVectorDbl("scratch/solution.mtx", X_in);

	assert(X.size() == X_in.size());

	Eigen::VectorXd diff = X-X_in;
	cout << "diff.norm() = " << diff.norm() << endl;

	printf("% 18.16e\n", (X.transpose()*H0.A*X)(0,0));
	printf("% 18.16e\n", (X_in.transpose()*H0.A*X_in)(0,0));
*/

//	for(int i=0; i<20; ++i)
//		printf("%5i: % f % f    -   % e\n", i, X[i], X_in[i], (X[i]-X_in[i])/X[i]);

}


/*
void test2(int N, dbl phi, int seed, bool write)
{
	//Create the system
	CStaticState<DIM> s(N);
	if(write)
	{
		s.RandomizePositions(seed);
		s.SetRadiiPolyUniform();
		s.SetPackingFraction(phi);

		dmat Trans;
		s.GetBox()->GetTransformation(Trans);
		Trans(1,0) += Trans(0,0)*0.1;
		s.GetBox()->SetTransformation(Trans);

		CStaticComputer<DIM> c(s);

		CStaticDatabase<DIM> db(N,"temp2.nc",NcFile::Replace);
		CStdDataDatabase<DIM> datadb("temp3.nc",NcFile::Replace);
		
		db.Write(s);

		CSimpleMinimizer<DIM> miner(c, CSimpleMinimizer<DIM>::FIRE);
		if(c.StdPrepareSystem())
		{
			printf("The bonds list is empty...\n");
			return;
		}
		c.CalculateStdData();
		c.Data.Print();
		
		db.Write(s, 4);
		datadb.Write(c.Data, 4);
	}else{
		CStaticDatabase<DIM> db(N,"temp2.nc",NcFile::ReadOnly);
		db.Read(s,4);
		
		CStaticComputer<DIM> c(s);
		CSimpleMinimizer<DIM> miner(c, CSimpleMinimizer<DIM>::FIRE);
		if(c.StdPrepareSystem())
		{
			printf("The bonds list is empty...\n");
			return;
		}
		c.CalculateStdData();
		c.Data.Print();
		
		CStdDataDatabase<DIM> datadb("temp3.nc",NcFile::ReadOnly);
		CStdData<DIM> data;
		datadb.Read(data, 4);
		data.Print();
	}
}
*/



/*


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


}
*/


int main(int argc, char* argv[])
{
	int N = 256;			//n
	dbl Lp = -1.0;			//p
	dbl phi = 0.845;			//f
	int r = 1;				//r
	int NFixedParticles = 0;//x
	bool write = false;		//w

	int c;
	while((c=getopt(argc, argv, "n:p:f:r:w")) != -1)
		switch(c)
		{
			case 'n':	N = atoi(optarg); break;
			case 'p':	Lp = atof(optarg); break;
			case 'f':	phi = atof(optarg); break;
			case 'r':	r = atoi(optarg); break;
			case 'w':	write = true; break;
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








