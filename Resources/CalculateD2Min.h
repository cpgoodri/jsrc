#ifndef CALC_D2MIN_H
#define CALC_D2MIN_H

#include "State/StaticState.h"
#include "Computers/StaticComputer.h"


namespace LiuJamming
{


template<int Dim>
class CD2MinCalculator
{
public:
	static void ComputeD2Min(CStaticComputer<Dim> &c, Eigen::VectorXd &dispVec, Eigen::VectorXd &D2Min);
	static void ComputeLocalStrain(CStaticComputer<Dim> &c, Eigen::VectorXd &dispVec, Eigen::VectorXd &LocalStrain);
};

template<int Dim>
void CD2MinCalculator<Dim>::ComputeD2Min(CStaticComputer<Dim> &c, Eigen::VectorXd &dispVec, Eigen::VectorXd &D2Min)
{
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;

	int NPp = c.Bonds.GetN();

	vector< vector<int> > bond_map(NPp, vector<int>());
	//for(typename vector< CBond<Dim> >::const_iterator b=c.Bonds.begin(); b!=c.Bonds.end(); ++b)
	for(int bi = 0; bi < c.Bonds.GetNBonds(); ++bi)
	{
		bond_map[c.Bonds[bi].i].push_back(bi);
		bond_map[c.Bonds[bi].j].push_back(bi);
	}
	assert(NPp == c.RattlerMap.size());
	assert(NPp*Dim == dispVec.size());

	//Initialize D2Min to zero
	D2Min = Eigen::VectorXd::Zero(NPp);

	dmat X, Y, Transform;
	dbl DispSqSum, tD2Min, ttD2Min;
	dvec dr0, dr1;
	dbl sgn;
	int bi, j;
	for(int i=0; i<NPp; ++i)
	{
		X = Y = Eigen::Matrix<dbl,Dim,Dim>::Zero();
		DispSqSum = 0.0;

		for(typename vector<int>::iterator bi_it = bond_map[i].begin(); bi_it != bond_map[i].end(); ++bi_it)
		{
			bi = (*bi_it);
			if(c.Bonds[bi].i == i)
			{
				j = c.Bonds[bi].j;
				sgn=1.;
			}
			else{
				assert(c.Bonds[bi].j == i);
				j = c.Bonds[bi].i;
				sgn=-1.;
			}

			dr0 = dr1 = sgn*c.Bonds[bi].r;
			for(int dd=0; dd<Dim; ++dd)
				dr1[dd] += sgn*(dispVec[Dim*j+dd]-dispVec[Dim*i+dd]);

			X += dr1*dr0.transpose();
			Y += dr0*dr0.transpose();
			DispSqSum += dr1.dot(dr1);
		}

		Transform = X*Y.inverse();
		//cout << Transform << endl;

		tD2Min = (DispSqSum - 2.0*(Transform*X.transpose()).trace() + (Transform * Y * Transform.transpose()).trace());
		D2Min[i] = tD2Min;
	}
}

template<int Dim>
void CD2MinCalculator<Dim>::ComputeLocalStrain(CStaticComputer<Dim> &c, Eigen::VectorXd &dispVec, Eigen::VectorXd &LocalStrain)
{
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;

	int NPp = c.Bonds.GetN();

	vector< vector<int> > bond_map(NPp, vector<int>());
	//for(typename vector< CBond<Dim> >::const_iterator b=c.Bonds.begin(); b!=c.Bonds.end(); ++b)
	for(int bi = 0; bi < c.Bonds.GetNBonds(); ++bi)
	{
		bond_map[c.Bonds[bi].i].push_back(bi);
		bond_map[c.Bonds[bi].j].push_back(bi);
	}
	assert(NPp == c.RattlerMap.size());
	assert(NPp*Dim == dispVec.size());

	//Initialize D2Min to zero
	LocalStrain = Eigen::VectorXd::Zero(NPp);

	dmat X, Y, Transform;
	dbl DispSqSum, tD2Min, ttD2Min;
	dvec dr0, dr1;
	dbl sgn;
	int bi, j;
	for(int i=0; i<NPp; ++i)
	{
		X = Y = Eigen::Matrix<dbl,Dim,Dim>::Zero();
		DispSqSum = 0.0;

		for(typename vector<int>::iterator bi_it = bond_map[i].begin(); bi_it != bond_map[i].end(); ++bi_it)
		{
			bi = (*bi_it);
			if(c.Bonds[bi].i == i)
			{
				j = c.Bonds[bi].j;
				sgn=1.;
			}
			else{
				assert(c.Bonds[bi].j == i);
				j = c.Bonds[bi].i;
				sgn=-1.;
			}

			dr0 = dr1 = sgn*c.Bonds[bi].r;
			for(int dd=0; dd<Dim; ++dd)
				dr1[dd] += sgn*(dispVec[Dim*j+dd]-dispVec[Dim*i+dd]);

			X += dr1*dr0.transpose();
			Y += dr0*dr0.transpose();
//			DispSqSum += dr1.dot(dr1);
		}

		Transform = X*Y.inverse();
		Transform = 0.5*(Transform - Transform.transpose());

		
		for(typename vector<int>::iterator bi_it = bond_map[i].begin(); bi_it != bond_map[i].end(); ++bi_it)
		{
			bi = (*bi_it);
			if(c.Bonds[bi].i == i)
			{
				j = c.Bonds[bi].j;
				sgn=1.;
			}
			else{
				assert(c.Bonds[bi].j == i);
				j = c.Bonds[bi].i;
				sgn=-1.;
			}

			dr0 = dr1 = sgn*c.Bonds[bi].r;
			for(int dd=0; dd<Dim; ++dd)
				dr1[dd] += sgn*(dispVec[Dim*j+dd]-dispVec[Dim*i+dd]);
			
			LocalStrain[i] = (dr1 - Transform*dr0).squaredNorm();
//			X += dr1*dr0.transpose();
//			Y += dr0*dr0.transpose();
//			DispSqSum += dr1.dot(dr1);
		}



		//tD2Min = (DispSqSum - 2.0*(Transform*X.transpose()).trace() + (Transform * Y * Transform.transpose()).trace());
		//D2Min[i] = tD2Min;
	}
}






}

#endif //CALC_D2MIN_H
