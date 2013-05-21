#ifndef STATICCOMPUTER

#define STATICCOMPUTER

/////////////////////////////////////////////////////////////////////////////////
//Static Computer class. 
//
//Description
//		This class maintains a spatial partition of particle systems so that
//		only local regions of a given particle need to be searched for neighbors.
//		The f 
//
//Variables
// 		Particle positions as an eigen vector
//		Particle radii as an eigen vector
//		Number of particles as an integer
//		Box as a CBox class
//		Potential as a CPotential class
//		
//Implements
//		Reading/Writing to netcdf files.
//		Getting/Setting particle positions.
//						particle radii.
//						box.
//						potential.
//
//
//File Formats
//		NetCDF
//			## Note: only files of a single N and Dimension may be stored in a given 
//			## file. To denote whether the NetCDF file has been populated with variables,
//			## dimensions, etc... an attribute "System_Populated" gets created.
//			-> Three dimensions: record number (unlimited), degrees of freedom (DOF), and
//			   particle number.
//			-> Dimension is stored as an attribute
//			-> Particle positions are stored as a variable.
//			-> Particle radii are stored as a variable.
//			-> Box is stored as a box (see CBox definition.)
//			-> Potential is stored as a potential (see CPotential definition.)
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include "../Potentials/Potentials.h"
#include "../Boundaries/Boxes.h"
#include "../Resources/MersenneTwister.h"
#include "../State/StaticState.h"
#include "BaseComputer.h"
#include "BondList.h"
#include "Grid.h"
#include "Data.h"
#include <list>


namespace LiuJamming
{

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
class CStaticComputer : public CBaseComputer<Dim>
{
private:
	typedef Eigen::Matrix<dbl,Dim,1>   dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	CStaticState<Dim> &State;
	CGrid<Dim> Grid;

public:
	CBondList<Dim> Bonds;	//!<Bond list for standard computations
	index_map RattlerMap;	//!<Map for rattlers.
	CStdData<Dim>  Data;	//!<storage for standard data

	CStaticComputer();
	CStaticComputer(CStaticState<Dim> &_State);
	CStaticComputer(const CStaticComputer<Dim> &_Copy);
	
		
//get and set the state
	const CStaticState<Dim> &GetState() const;
	void SetState(CStaticState<Dim> & State);

//Get generic info from State
	dbl GetVolume() const;
	int GetNdof() const;

//Compute the bond list
	void ComputeBondList(CBondList<Dim> &bonds);
	void ComputeBondList_NoGrid(CBondList<Dim> &bonds) const;

	void StdPrepareSystem();
	void CalculateStdData(bool CalcCijkl = true, bool CalcHess = true);
	void CalculateStdData(CStdData<Dim> &data, bool CalcCijkl=true, bool CalcHess = true);
	void CalculateStdData_Unstressed(bool CalcCijkl=true, bool CalcHess = true);
	void CalculateStdData_Unstressed(CStdData<Dim> &data, bool CalcCijkl=true, bool CalcHess = true);

	void ComputeHessian(Eigen::SparseMatrix<dbl> &hess) const;
	void ComputeHessian(CBondList<Dim> const &bonds, Eigen::SparseMatrix<dbl> &hess) const;
	void ComputeHessianBZ(Eigen::SparseMatrix<cdbl> &hess, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const;
	void ComputeHessianBZ(CBondList<Dim> const &bonds, Eigen::SparseMatrix<cdbl> &hess, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const;

/*
//Compute the energy of the system
    dbl ComputeEnergy();
    
//compute the gradient of the energy of the system.
    void ComputeGradient(Eigen::VectorXd &tar);
    void ComputeForce(Eigen::VectorXd &tar);
    
//Hessian stuff
	void ComputeHessian(Eigen::MatrixXd &tar);

//Dynamical Matrix Stuff (all done in mass-normalized coordinates)
//computes at q = 0
	void ComputeDynamicalMatrix(Eigen::MatrixXd &tar);
*/

//Needed for minimization routines
	void Evaluate(Eigen::VectorXd &grad, dbl &fx);
	bool Progress(Eigen::VectorXd const &grad, dbl fx, int iteration, dbl tol);
	void Move(Eigen::VectorXd const &step);
	dbl GetMinimizationTimeScale() const;
    

	void CalculateBZTransformation(Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const;
	void CalculateBZTransformation(CBondList<Dim> const &bonds, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const;
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION    //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CStaticComputer<Dim>::CStaticComputer()
{
}


template <int Dim>
CStaticComputer<Dim>::CStaticComputer(CStaticState<Dim> &_State) : State(_State), Grid(&State)
{	
	Grid.Allocate();
}

template <int Dim>
CStaticComputer<Dim>::CStaticComputer(const CStaticComputer<Dim> &_Copy) : State(_Copy.State), Grid(&State)
{
	Grid.Allocate();
}


//get and set the state
template <int Dim>
const CStaticState<Dim> &CStaticComputer<Dim>::GetState() const
{
	return State;
}

template <int Dim>
void CStaticComputer<Dim>::SetState(CStaticState<Dim> & _State)
{
	State = _State;
	Grid.SetState(&State);
}
	
template <int Dim>
inline dbl CStaticComputer<Dim>::GetVolume() const
{
	return State.GetVolume();
}

template <int Dim>
inline int CStaticComputer<Dim>::GetNdof() const
{
	return State.GetParticleNumber()*Dim;
}

template <int Dim>
inline dbl CStaticComputer<Dim>::GetMinimizationTimeScale() const
{
	dbl avgRad = State.GetAverageRadius();
	return avgRad; //This is emperical!!! NOT GENERAL!
}

//Compute the bond list
template <int Dim>
void CStaticComputer<Dim>::ComputeBondList(CBondList<Dim> &bonds)
{
	Grid.Construct();
	//Grid.PrintGrid();

	bonds.SetN(State.GetParticleNumber());
	dvec Displacement;
	dbl sigma, rlen, rlen2, E, g, k;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++) //(*it) is the particle index of a potential neighbor
		{
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				sigma = State.GetRadius(i) + State.GetRadius((*it));
				rlen2 = Displacement.squaredNorm();
				if(rlen2 < sigma*sigma)
				{
					rlen = sqrt(rlen2);
					State.GetPotential()->ComputeDerivatives012(rlen, sigma, E, g, k);
					bonds.AddBond( CBond<Dim>(i, (*it), rlen, E, g, k, Displacement) );
					//bonds.AddBond( CBond<Dim>(i, (*it), sigma, rlen, E, g, k, Displacement) );
				}
			}
		}
	}
	bonds.SetVolume(GetVolume());
	
}

template <int Dim>
void CStaticComputer<Dim>::ComputeBondList_NoGrid(CBondList<Dim> &bonds) const
{
	dvec Displacement;
	bonds.SetN(State.GetParticleNumber());
	dbl sigma, rlen, rlen2, E, g, k;
	for(int i = 0 ; i < State.GetParticleNumber()-1 ; i++)
		for(int j = i+1 ; j < State.GetParticleNumber() ; j++)
		{
			State.GetDisplacement(i,j,Displacement);
			sigma = State.GetRadius(i) + State.GetRadius(j);
			rlen2 = Displacement.squaredNorm();
			if(rlen2 < sigma*sigma)
			{
				rlen = sqrt(rlen2);
				State.GetPotential()->ComputeDerivatives012(rlen, sigma, E, g, k);
				bonds.AddBond( CBond<Dim>(i, j, rlen, E, g, k, Displacement) );
				//bonds.AddBond( CBond<Dim>(i, j, sigma, rlen, E, g, k, Displacement) );
			}
		}
	//bonds.Volume = 0.;
	bonds.SetVolume(GetVolume());
}










template <int Dim>
void CStaticComputer<Dim>::StdPrepareSystem()
{
	ComputeBondList(Bonds);				//Use the member variable derived from CBaseComputer
	Bonds.RemoveRattlers(RattlerMap);	//Remove rattlers
}

template <int Dim>
void CStaticComputer<Dim>::CalculateStdData(bool CalcCijkl, bool CalcHess)
{
	CalculateStdData(Data,CalcCijkl,CalcHess);
}
template <int Dim>
void CStaticComputer<Dim>::CalculateStdData(CStdData<Dim> &data, bool CalcCijkl, bool CalcHess)
{
	Bonds.CalculateStdData(data, CalcCijkl,CalcHess);
}


template <int Dim>
void CStaticComputer<Dim>::CalculateStdData_Unstressed(bool CalcCijkl, bool CalcHess)
{
	CalculateStdData_Unstressed(Data,CalcCijkl,CalcHess);
}


template <int Dim>
void CStaticComputer<Dim>::CalculateStdData_Unstressed(CStdData<Dim> &data, bool CalcCijkl, bool CalcHess)
{
	CBondList<Dim> BondsTemp = Bonds;
	BondsTemp.MakeUnstressed();
	BondsTemp.CalculateStdData(data,CalcCijkl,CalcHess);
}











//Needed for minimization routines
template <int Dim>
void CStaticComputer<Dim>::Evaluate(Eigen::VectorXd &grad, dbl &fx) 
{
	CBondList<Dim> bonds;
	ComputeBondList(bonds);
	fx = bonds.ComputeGradient(grad);
};

template <int Dim>
bool CStaticComputer<Dim>::Progress(Eigen::VectorXd const &grad, dbl fx, int iteration, dbl tol) 
{
	dbl gradNorm = grad.norm();
	dbl max_grad = max_abs_element(grad.size(), grad.data());
	bool converged = (max_grad < tol)?true:false;
	if(converged){ printf("converged\n"); fflush(stdout);}
	if(iteration%1000==0 || converged)
	{
		printf("% 11i   %22.20e   % e   % e\n", iteration, fx, gradNorm, max_grad);
		fflush(stdout);
	}
	return converged;
};

template <int Dim>
inline void CStaticComputer<Dim>::Move(Eigen::VectorXd const &step) 
{
	State.MoveParticles(step);
};
  








template <int Dim>
void CStaticComputer<Dim>::ComputeHessian(Eigen::SparseMatrix<dbl> &hess) const
{
	ComputeHessian(Bonds, hess);
}

template <int Dim>
void CStaticComputer<Dim>::ComputeHessian(CBondList<Dim> const &bonds, Eigen::SparseMatrix<dbl> &hess) const
{
	bonds.ComputeHessian(hess);
}

template <int Dim>
void CStaticComputer<Dim>::ComputeHessianBZ(Eigen::SparseMatrix<cdbl> &hess, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const
{
	ComputeHessianBZ(Bonds, hess, transformation, k);
}
template <int Dim>
void CStaticComputer<Dim>::ComputeHessianBZ(CBondList<Dim> const &bonds, Eigen::SparseMatrix<cdbl> &hess, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const
{
	CalculateBZTransformation(bonds, transformation, k);
	Eigen::SparseMatrix<dbl> H0;
	ComputeHessian(bonds, H0);
	printf("transformation -> %i %i\n", transformation.rows(), transformation.cols());
	printf("H0             -> %i %i\n", H0.rows(), H0.cols());
	hess = transformation.adjoint() * H0 * transformation;
}

template <int Dim>
void CStaticComputer<Dim>::CalculateBZTransformation(Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const
{
	CalculateBZTransformation(Bonds, transformation, k);
}

template <int Dim>
void CStaticComputer<Dim>::CalculateBZTransformation(CBondList<Dim> const &bonds, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const
{
	//Get the real positions
	Eigen::VectorXd RealPos;
	State.GetPositions(RealPos);

	//Remove Positions corresponding to rattlers
	Eigen::VectorXd RealPosNoRatt = Eigen::VectorXd::Zero(Dim*RattlerMap.size());
	for(int im=0; im<RattlerMap.size(); ++im)
		RealPosNoRatt.segment<Dim>(Dim*im) = RealPos.segment<Dim>(Dim*RattlerMap[im]);
	bonds.PlaneWaveAnsatz(transformation, RealPosNoRatt, k);
}




/*
//Compute the energy of the system
template <int Dim>
dbl CStaticComputer<Dim>::ComputeEnergy()
{
	Grid.Construct();

	Eigen::Matrix<dbl,Dim,1> Displacement;
	dbl Energy = 0.0;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++)
		{
			cout.flush();
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				dbl rij = State.GetRadius(i) + State.GetRadius((*it));
				dbl sqn = Displacement.squaredNorm();
				if(sqn<rij*rij)
					Energy+= State.GetPotential()->Compute(sqrt(sqn),rij);
			}
		}
	}


	return Energy;
}
    
//compute the gradient of the energy of the system.
template <int Dim>
void CStaticComputer<Dim>::ComputeGradient(Eigen::VectorXd &tar)
{
	Grid.Construct();

	Eigen::Matrix<dbl,Dim,1> Displacement;
	dbl Energy = 0.0;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++)
		{
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				dbl rij = State.GetRadius(i) + State.GetRadius((*it));
				dbl sqn = Displacement.squaredNorm();
				if(Displacement.squaredNorm()<rij*rij){
					sqn = sqrt(sqn);
					dbl dU = State.GetPotential()->ComputeFirstDerivative(sqn,rij);
					for(int j = 0; j<Dim ; j++)
					{
						tar(Dim*i+j) += -dU*Displacement(j)/rij;
                        tar(Dim*(*it)+j) += dU*Displacement(j)/rij;
					}
				}
			}
		}
	}
}

template <int Dim>
void CStaticComputer<Dim>::ComputeForce(Eigen::VectorXd &tar)
{
	ComputeGradient(tar);
	tar = -tar;
}

    
//Hessian stuff
template <int Dim>
void CStaticComputer<Dim>::ComputeHessian(Eigen::MatrixXd &tar) 
{
	Grid.Construct();

	Eigen::Matrix<dbl,Dim,1> Displacement;
	dbl Energy = 0.0;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++)
		{
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				dbl rij = State.GetRadius(i) + State.GetRadius((*it));
				dbl sqn = Displacement.squaredNorm();
				if(Displacement.squaredNorm()<rij*rij){
					sqn = sqrt(sqn);
					dbl dU = State.GetPotential()->ComputeFirstDerivative(sqn,rij);					
					dbl d2U = State.GetPotential()->ComputeSecondDerivative(sqn,rij);
					//map the block hessian into the actual hessian
					int j = (*it);
                    for(int q = 0 ; q<Dim ; q++)
					{
						dbl ii = dU*(1.0-Displacement(q)*Displacement(q)/sqn)/rij + d2U*Displacement(q)*Displacement(q)/sqn;
						//i_qi_q
						tar(Dim*i+q,Dim*i+q)+=ii;
						tar(Dim*j+q,Dim*j+q)+=ii;
						//i_qj_q
						tar(Dim*i+q,Dim*j+q)-=ii;
						tar(Dim*j+q,Dim*i+q)-=ii;
						for(int q2 = q+1 ; q2<Dim ; q2++)
						{
							dbl ij = Displacement(q)*Displacement(q2)/sqn*(-dU/rij + d2U);
						
							//i_qi_q2
							tar(Dim*i+q,Dim*i+q2)+=ij;
							tar(Dim*i+q2,Dim*i+q)+=ij;
							tar(Dim*j+q,Dim*j+q2)+=ij;
							tar(Dim*j+q2,Dim*j+q)+=ij;
							
							//i_qj_q2
							tar(Dim*i+q,Dim*j+q2)-=ij;
							tar(Dim*i+q2,Dim*j+q)-=ij;
							tar(Dim*j+q,Dim*i+q2)-=ij;
							tar(Dim*j+q2,Dim*i+q)-=ij;
						}
					}
				}
			}
		}
	}
}

//Dynamical Matrix Stuff (all done in mass-normalized coordinates)
//computes at q = 0
template <int Dim>
void CStaticComputer<Dim>::ComputeDynamicalMatrix(Eigen::MatrixXd &tar)
{
	Grid.Construct();

	Eigen::Matrix<dbl,Dim,1> Displacement;
	dbl Energy = 0.0;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++)
		{
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				dbl rij = State.GetRadius(i) + State.GetRadius((*it));
				dbl sqn = Displacement.squaredNorm();
				if(Displacement.squaredNorm()<rij*rij){
					dbl mass1 = 1.0/pow(State.GetRadius(i)*State.GetRadius(i),Dim/2.0);//1.0/Particles[i].Radius/Particles[i].Radius/Scale/Scale;
                    dbl mass2 = 1.0/pow(State.GetRadius((*it))*State.GetRadius((*it)),Dim/2.0);//1.0/Particles[j].Radius/Particles[j].Radius/Scale/Scale;
                    dbl mass12 = 1.0/pow(State.GetRadius(i)*State.GetRadius((*it)),Dim/2.0);
					sqn = sqrt(sqn);
					dbl dU = State.GetPotential()->ComputeFirstDerivative(sqn,rij);					
					dbl d2U = State.GetPotential()->ComputeSecondDerivative(sqn,rij);
					//map the block hessian into the actual hessian
					int j = (*it);
                    for(int q = 0 ; q<Dim ; q++)
					{
						dbl ii = dU*(1.0-Displacement(q)*Displacement(q)/sqn)/rij + d2U*Displacement(q)*Displacement(q)/sqn;
						//i_qi_q
						tar(Dim*i+q,Dim*i+q)+=mass1*ii;
						tar(Dim*j+q,Dim*j+q)+=mass2*ii;
						//i_qj_q
						tar(Dim*i+q,Dim*j+q)-=mass12*ii;
						tar(Dim*j+q,Dim*i+q)-=mass12*ii;
						for(int q2 = q+1 ; q2<Dim ; q2++)
						{
							dbl ij = Displacement(q)*Displacement(q2)/sqn*(-dU/rij + d2U);
						
							//i_qi_q2
							tar(Dim*i+q,Dim*i+q2)+=mass1*ij;
							tar(Dim*i+q2,Dim*i+q)+=mass1*ij;
							tar(Dim*j+q,Dim*j+q2)+=mass2*ij;
							tar(Dim*j+q2,Dim*j+q)+=mass2*ij;
							
							//i_qj_q2
							tar(Dim*i+q,Dim*j+q2)-=mass12*ij;
							tar(Dim*i+q2,Dim*j+q)-=mass12*ij;
							tar(Dim*j+q,Dim*i+q2)-=mass12*ij;
							tar(Dim*j+q2,Dim*i+q)-=mass12*ij;
						}
					}
				}
			}
		}
	}
}
    */

}

#endif
