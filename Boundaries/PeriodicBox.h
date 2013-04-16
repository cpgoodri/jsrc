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

template <int Dim, int NonPeriodicDim=0>
class CPeriodicBox : public CBox<Dim>
{	
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;

public:
//constructors and copy operators
	CPeriodicBox();
	CPeriodicBox(const dmat Trans);
	CPeriodicBox(dbl sx, dbl sy, dbl sz);
	CPeriodicBox(const CPeriodicBox &box);
	
	const CPeriodicBox<Dim,NonPeriodicDim> &operator=(const CPeriodicBox<Dim,NonPeriodicDim> &box);

//functions to write box configurations
	string DataToString();
	
//functions to read box configurations
 	void StringToData(string Data);
 	CBox<Dim> *Create();
 	
//functions involving the boundary
	void MoveParticles(Eigen::VectorXd &Points, const Eigen::VectorXd &Displacements);
	void ApplyPeriodicBC(Eigen::VectorXd &Points);
	void MinimumDisplacement(const dvec &PointA, const dvec &PointB, dvec &Displacement) const;
	void MinimumDisplacement(const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointA, const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointB, dvec &Displacement) const;
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//constructors and copy operators
template <int Dim, int NonPeriodicDim>
CPeriodicBox<Dim,NonPeriodicDim>::CPeriodicBox() : CBox<Dim>()
{

}

template <int Dim, int NonPeriodicDim>
CPeriodicBox<Dim,NonPeriodicDim>::CPeriodicBox(const dmat Trans) : CBox<Dim>(Trans)
{


}

template <int Dim, int NonPeriodicDim>
CPeriodicBox<Dim,NonPeriodicDim>::CPeriodicBox(dbl sx, dbl sy, dbl sz) : CBox<Dim>(sx,sy,sz)
{

}

template <int Dim, int NonPeriodicDim>
CPeriodicBox<Dim,NonPeriodicDim>::CPeriodicBox(const CPeriodicBox &box) : CBox<Dim>(box)
{

}

template <int Dim, int NonPeriodicDim>	
const CPeriodicBox<Dim,NonPeriodicDim> &CPeriodicBox<Dim,NonPeriodicDim>::operator=(const CPeriodicBox<Dim,NonPeriodicDim> &box)
{
	CBox<Dim>::operator=(box);
	return *this;
}

//functions to write box configurations
template <int Dim, int NonPeriodicDim>
string CPeriodicBox<Dim,NonPeriodicDim>::DataToString()
{
	return "PeriodicBox";
}
	
//functions to read box configurations
template <int Dim, int NonPeriodicDim>
void CPeriodicBox<Dim,NonPeriodicDim>::StringToData(string Data)
{

}

template <int Dim, int NonPeriodicDim>
CBox<Dim> *CPeriodicBox<Dim,NonPeriodicDim>::Create()
{
	return new CPeriodicBox<Dim>();
}


//	class PeriodicBCs
//		This class implements a function to move particle positions
//			according to periodic boundary conditions
//		It is a templated class with a partial specialization
//			to handle the common case of full periodic BCs.
template<int Dim, int NonPeriodicPos>class PeriodicBCs
{
	static inline void apply(Eigen::VectorXd &Points){
		assert(Points.cols()%Dim==0);
		int np = Points.cols()/Dim;
		for(int i=0; i<np; i++)
			for(int dd=NonPeriodicDim; dd<Dim; dd++)
				Points(i*Dim+dd) -= floor(Points(i*Dim+dd));
	};
};
template<int Dim>class PeriodicBCs<Dim,0>
{
	static inline void apply(Eigen::VectorXd &Points){
		for(int i = 0; i< Points.cols() ; i++)
			Points(i) -= floor(Points(i));
	};
};

//functions involving the boundary
template <int Dim, int NonPeriodicDim>
void CPeriodicBox<Dim,NonPeriodicDim>::MoveParticles(Eigen::VectorXd &Points,const Eigen::VectorXd &Displacements)
{
	Points += Displacements;
	PeriodicBCs<Dim,NonPeriodicPos>::apply(Points);
//	ApplyPeriodicBC<Dim,NonPeriodicDim>(Points);
}

//These two methods assume:
//	1. the box has dimensions 1.
//	2. all particles are within the box.
//
//I'm not sure if this will work with free boundary conditions
template <int Dim, int NonPeriodicDim>
void CPeriodicBox<Dim,NonPeriodicDim>::MinimumDisplacement(const dvec &PointA,const dvec &PointB, dvec &Displacement) const
{
	Displacement = PointA-PointB;
	for(int i = NonPeriodicDim ; i < Dim ; i++)
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
}

template <int Dim, int NonPeriodicDim>
void CPeriodicBox<Dim,NonPeriodicDim>::MinimumDisplacement(const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointA,const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointB, dvec &Displacement) const
{
	Displacement = PointA-PointB;
	for(int i = NonPeriodicDim ; i < Dim ; i++)
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
}

}

#endif
