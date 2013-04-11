#ifndef PERIODICBOX

#define PERIODICBOX

/////////////////////////////////////////////////////////////////////////////////
//Periodic Box class. 
//
//Description
//		This is a class that defines periodic boundary conditions. It inherits from
//		the CBox class. The plaintext name is "PeriodicBox".
//
//Variables
// 		
//		
//Implements
//		Computes minimal distances using periodic boundary conditions.
//		Displaces particles while respecting the periodic boundary conditions.
//		Maps the box data to a string.
//
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include "../Resources/MersenneTwister.h"
#include "Box.h"
#include <map>
#include <string>

namespace LiuJamming
{

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
class CPeriodicBox : public CBox<Dim>
{	
public:
//constructors and copy operators
	CPeriodicBox();
	CPeriodicBox(const Eigen::Matrix<double,Dim,Dim> Trans);
	CPeriodicBox(double sx, double sy, double sz);
	CPeriodicBox(const CPeriodicBox &box);
	
	const CPeriodicBox<Dim> &operator=(const CPeriodicBox<Dim> &box);

//functions to write box configurations
	string DataToString();
	
//functions to read box configurations
 	void StringToData(string Data);
 	CBox<Dim> *Create();
 	
//functions involving the boundary
	void MoveParticles(Eigen::VectorXd &Points,const Eigen::VectorXd &Displacements);
	void MinimumDisplacement(const Eigen::Matrix<double,Dim,1> &PointA,const Eigen::Matrix<double,Dim,1> &PointB, Eigen::Matrix<double,Dim,1> &Displacement) const;
	void MinimumDisplacement(const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointA,const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointB, Eigen::Matrix<double,Dim,1> &Displacement) const;
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//constructors and copy operators
template <int Dim>
CPeriodicBox<Dim>::CPeriodicBox() : CBox<Dim>()
{

}

template <int Dim>
CPeriodicBox<Dim>::CPeriodicBox(const Eigen::Matrix<double,Dim,Dim> Trans) : CBox<Dim>(Trans)
{


}

template <int Dim>
CPeriodicBox<Dim>::CPeriodicBox(double sx, double sy, double sz) : CBox<Dim>(sx,sy,sz)
{

}

template <int Dim>
CPeriodicBox<Dim>::CPeriodicBox(const CPeriodicBox &box) : CBox<Dim>(box)
{

}

template <int Dim>	
const CPeriodicBox<Dim> &CPeriodicBox<Dim>::operator=(const CPeriodicBox<Dim> &box)
{
	CBox<Dim>::operator=(box);
	return *this;
}

//functions to write box configurations
template <int Dim>
string CPeriodicBox<Dim>::DataToString()
{
	return "PeriodicBox";
}
	
//functions to read box configurations
template <int Dim>
void CPeriodicBox<Dim>::StringToData(string Data)
{

}

template <int Dim>
CBox<Dim> *CPeriodicBox<Dim>::Create()
{
	return new CPeriodicBox<Dim>();
}
 	
//functions involving the boundary
template <int Dim>
void CPeriodicBox<Dim>::MoveParticles(Eigen::VectorXd &Points,const Eigen::VectorXd &Displacements)
{
	Points += Displacements;
	for(int i = 0; i< Points.cols() ; i++)
		Points(i) -= floor(Points(i));
}

template <int Dim>
void CPeriodicBox<Dim>::MinimumDisplacement(const Eigen::Matrix<double,Dim,1> &PointA,const Eigen::Matrix<double,Dim,1> &PointB, Eigen::Matrix<double,Dim,1> &Displacement) const
{
	Displacement = PointA-PointB;
	for(int i = 0 ; i < Dim ; i++)
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
}

template <int Dim>
void CPeriodicBox<Dim>::MinimumDisplacement(const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointA,const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointB, Eigen::Matrix<double,Dim,1> &Displacement) const
{
	Displacement = PointA-PointB;
	for(int i = 0 ; i < Dim ; i++)
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
}

}

#endif
