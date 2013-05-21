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
#include "BaseBox.h"
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
	typedef Eigen::VectorBlock<Eigen::VectorXd,Dim> dvecBlock;

public:
//constructors and copy operators
	CPeriodicBox();
	CPeriodicBox(const dmat Trans);
	CPeriodicBox(dbl sx, dbl sy, dbl sz);
	CPeriodicBox(const CPeriodicBox &box);
	
	const CPeriodicBox<Dim,NonPeriodicDim> &operator=(const CPeriodicBox<Dim,NonPeriodicDim> &box);

//functions to read and write box configurations
	static string GetName();
	string DataToString() const;
 	void StringToData(string Data);
 	CBox<Dim> *Create();
	
//get a list of the periodic dimensions
	void GetPeriodicDimensions(std::vector<int> &) const;
 	
//functions involving the boundary
	void MoveParticles(Eigen::VectorXd &Points, const Eigen::VectorXd &Displacements) const;
	void MoveParticle(dvecBlock Point, dvec const &Displacement) const;
	void MinimumDisplacement(const dvec &PointA, const dvec &PointB, dvec &Displacement) const;
	void MinimumDisplacement(const dvecBlock &PointA, const dvecBlock &PointB, dvec &Displacement) const;
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

//functions to read and write box configurations
template <int Dim, int NonPeriodicDim>
string CPeriodicBox<Dim,NonPeriodicDim>::GetName()
{
	stringstream name;
	name << "PeriodicBox_" << Dim << "D_" << NonPeriodicDim << "NonPDim";
	return name.str();
}
	
template <int Dim, int NonPeriodicDim>
string CPeriodicBox<Dim,NonPeriodicDim>::DataToString() const
{
	stringstream ss;
	ss << GetName() << ":" << CBox<Dim>::BoxSymmetry;
	return ss.str();
}
	
template <int Dim, int NonPeriodicDim>
void CPeriodicBox<Dim,NonPeriodicDim>::StringToData(string Data)
{
	vector<string> split = SplitString(Data,":");
	CBox<Dim>::BoxSymmetry = atoi(split[1].c_str());
}

template <int Dim, int NonPeriodicDim>
CBox<Dim> *CPeriodicBox<Dim,NonPeriodicDim>::Create()
{
	return new CPeriodicBox<Dim,NonPeriodicDim>();
}


//	class PeriodicBCs
//		This class implements a function to move particle positions
//			according to periodic boundary conditions
//		It is a templated class with a partial specialization
//			to handle the common case of full periodic BCs.
template<int Dim, int NonPeriodicDim> class PeriodicBCs
{
public:
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::VectorBlock<Eigen::VectorXd,Dim> dvecBlock;
	static inline void apply_particle(dvecBlock Point)
	{
		for(int dd=NonPeriodicDim; dd<Dim; dd++)
			Point(dd) -= floor(Point(dd));
	}
	static inline void apply_vector(Eigen::VectorXd &Points){
		assert(Points.cols()%Dim==0);
		int np = Points.cols()/Dim;
		for(int i=0; i<np; i++)
			for(int dd=NonPeriodicDim; dd<Dim; dd++)
				Points(i*Dim+dd) -= floor(Points(i*Dim+dd));
	};
};
template<int Dim>class PeriodicBCs<Dim,0>
{
public:
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::VectorBlock<Eigen::VectorXd,Dim> dvecBlock;
	static inline void apply_particle(dvecBlock Point)
	{
		for(int dd=0; dd<Dim; dd++)
			Point(dd) -= floor(Point(dd));
	}
	static inline void apply_vector(Eigen::VectorXd &Points){
		for(int i = 0; i< Points.cols() ; i++)
			Points(i) -= floor(Points(i));
	};
};

//functions involving the boundary
template <int Dim, int NonPeriodicDim>
inline void CPeriodicBox<Dim,NonPeriodicDim>::MoveParticles(Eigen::VectorXd &Points,const Eigen::VectorXd &Displacements) const
{
	Points += Displacements;
	PeriodicBCs<Dim,NonPeriodicDim>::apply_vector(Points);
}

template <int Dim, int NonPeriodicDim>
inline void CPeriodicBox<Dim,NonPeriodicDim>::MoveParticle(dvecBlock Point, dvec const &Displacement) const
{
	Point += Displacement;
	PeriodicBCs<Dim,NonPeriodicDim>::apply_particle(Point);
}




//These two methods assume:
//	1. the box has dimensions 1.
//	2. all particles are within the box.
//
//I'm not sure if this will work with free boundary conditions
template <int Dim, int NonPeriodicDim>
inline void CPeriodicBox<Dim,NonPeriodicDim>::MinimumDisplacement(const dvec &PointA,const dvec &PointB, dvec &Displacement) const
{
	Displacement = PointA-PointB;
	for(int i = NonPeriodicDim ; i < Dim ; i++)
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
}

template <int Dim, int NonPeriodicDim>
inline void CPeriodicBox<Dim,NonPeriodicDim>::MinimumDisplacement(const dvecBlock &PointA,const dvecBlock &PointB, dvec &Displacement) const
{
	Displacement = PointA-PointB;
	for(int i = NonPeriodicDim ; i < Dim ; i++)
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
}

template <int Dim, int NonPeriodicDim>
void CPeriodicBox<Dim,NonPeriodicDim>::GetPeriodicDimensions(std::vector<int> &pdims) const
{
	pdims.clear();
	for(int i=NonPeriodicDim; i<Dim; ++i)
		pdims.push_back(i);
}

}

#endif
