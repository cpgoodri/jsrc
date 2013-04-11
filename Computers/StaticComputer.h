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
#include "../Potentials/Potential.h"
#include "../Boundaries/Box.h"
#include "../Resources/MersenneTwister.h"
#include "../State/StaticState.h"
#include "BaseComputer.h"
#include "Grid.h"
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
	CStaticState<Dim> &State;
	CGrid<Dim> Grid;

public:
	CStaticComputer();
	CStaticComputer(CStaticState<Dim> &_State);
	CStaticComputer(const CStaticComputer<Dim> &_Copy);
	
		
//get and set the state
	const CStaticState<Dim> &GetState() const;
	void SetState(CStaticState<Dim> & State);
	

//Compute the energy of the system
    double ComputeEnergy();
    
//compute the gradient of the energy of the system.
    void ComputeGradient(Eigen::VectorXd &tar);
    void ComputeForce(Eigen::VectorXd &tar);
    
//Hessian stuff
    void ComputeHessian(Eigen::MatrixXd &tar);

//Dynamical Matrix Stuff (all done in mass-normalized coordinates)
//computes at q = 0
    void ComputeDynamicalMatrix(Eigen::MatrixXd &tar);
    
    
};

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
	

//Compute the energy of the system
template <int Dim>
double CStaticComputer<Dim>::ComputeEnergy()
{
	Grid.Construct();

	Eigen::Matrix<double,Dim,1> Displacement;
	double Energy = 0.0;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++)
		{
			cout.flush();
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				double rij = State.GetRadius(i) + State.GetRadius((*it));
				double sqn = Displacement.squaredNorm();
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

	Eigen::Matrix<double,Dim,1> Displacement;
	double Energy = 0.0;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++)
		{
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				double rij = State.GetRadius(i) + State.GetRadius((*it));
				double sqn = Displacement.squaredNorm();
				if(Displacement.squaredNorm()<rij*rij){
					sqn = sqrt(sqn);
					double dU = State.GetPotential()->ComputeFirstDerivative(sqn,rij);
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

	Eigen::Matrix<double,Dim,1> Displacement;
	double Energy = 0.0;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++)
		{
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				double rij = State.GetRadius(i) + State.GetRadius((*it));
				double sqn = Displacement.squaredNorm();
				if(Displacement.squaredNorm()<rij*rij){
					sqn = sqrt(sqn);
					double dU = State.GetPotential()->ComputeFirstDerivative(sqn,rij);					
					double d2U = State.GetPotential()->ComputeSecondDerivative(sqn,rij);
					//map the block hessian into the actual hessian
					int j = (*it);
                    for(int q = 0 ; q<Dim ; q++)
					{
						double ii = dU*(1.0-Displacement(q)*Displacement(q)/sqn)/rij + d2U*Displacement(q)*Displacement(q)/sqn;
						//i_qi_q
						tar(Dim*i+q,Dim*i+q)+=ii;
						tar(Dim*j+q,Dim*j+q)+=ii;
						//i_qj_q
						tar(Dim*i+q,Dim*j+q)-=ii;
						tar(Dim*j+q,Dim*i+q)-=ii;
						for(int q2 = q+1 ; q2<Dim ; q2++)
						{
							double ij = Displacement(q)*Displacement(q2)/sqn*(-dU/rij + d2U);
						
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

	Eigen::Matrix<double,Dim,1> Displacement;
	double Energy = 0.0;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++)
		{
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				double rij = State.GetRadius(i) + State.GetRadius((*it));
				double sqn = Displacement.squaredNorm();
				if(Displacement.squaredNorm()<rij*rij){
					double mass1 = 1.0/pow(State.GetRadius(i)*State.GetRadius(i),Dim/2.0);//1.0/Particles[i].Radius/Particles[i].Radius/Scale/Scale;
                    double mass2 = 1.0/pow(State.GetRadius((*it))*State.GetRadius((*it)),Dim/2.0);//1.0/Particles[j].Radius/Particles[j].Radius/Scale/Scale;
                    double mass12 = 1.0/pow(State.GetRadius(i)*State.GetRadius((*it)),Dim/2.0);
					sqn = sqrt(sqn);
					double dU = State.GetPotential()->ComputeFirstDerivative(sqn,rij);					
					double d2U = State.GetPotential()->ComputeSecondDerivative(sqn,rij);
					//map the block hessian into the actual hessian
					int j = (*it);
                    for(int q = 0 ; q<Dim ; q++)
					{
						double ii = dU*(1.0-Displacement(q)*Displacement(q)/sqn)/rij + d2U*Displacement(q)*Displacement(q)/sqn;
						//i_qi_q
						tar(Dim*i+q,Dim*i+q)+=mass1*ii;
						tar(Dim*j+q,Dim*j+q)+=mass2*ii;
						//i_qj_q
						tar(Dim*i+q,Dim*j+q)-=mass12*ii;
						tar(Dim*j+q,Dim*i+q)-=mass12*ii;
						for(int q2 = q+1 ; q2<Dim ; q2++)
						{
							double ij = Displacement(q)*Displacement(q2)/sqn*(-dU/rij + d2U);
						
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

    

}

#endif
