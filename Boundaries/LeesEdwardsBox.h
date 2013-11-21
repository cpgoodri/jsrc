#ifndef LEES_EDWARDS_BOX
#define LEES_EDWARDS_BOX

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

//!Class to implement a periodic box
template <int Dim>
class CLeesEdwardsBox : public CBox<Dim>
{	
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	typedef Eigen::VectorBlock<Eigen::VectorXd,Dim> dvecBlock;

	int ShearAxis;
	int ShearDirection;
	dbl Strain;

public:
//constructors and copy operators
	CLeesEdwardsBox(int sa=1, int sd=0, dbl s=0.);
	CLeesEdwardsBox(const dmat &Trans, int sa=1, int sd=0, dbl s=0.);
	CLeesEdwardsBox(const CLeesEdwardsBox &box);
	
	const CLeesEdwardsBox<Dim> &operator=(const CLeesEdwardsBox<Dim> &box);

//functions to read and write box configurations
	static string GetName();
	virtual string DataToString() const;
	virtual void StringToData(string Data);
	virtual CLeesEdwardsBox<Dim>* Clone() const;
	
//get a list of the periodic dimensions
	void GetPeriodicDimensions(std::vector<int> &) const;
 	
//functions involving the boundary
	virtual void MoveParticles(Eigen::VectorXd &Points, const Eigen::VectorXd &Displacements) const;
	virtual void MoveParticle(dvecBlock Point, dvec const &Displacement) const;
	virtual void MinimumDisplacement(const dvec &PointA, const dvec &PointB, dvec &Displacement) const;
	virtual void MinimumDisplacement(const dvecBlock &PointA, const dvecBlock &PointB, dvec &Displacement) const;
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
CLeesEdwardsBox<Dim>::CLeesEdwardsBox(int sa=1, int sd=0, dbl s=0.) 
	: CBox<Dim>(), StrainAxis(sa), StrainDirection(sd), Strain(s)
{}

template <int Dim>
CLeesEdwardsBox<Dim>::CLeesEdwardsBox(const dmat &Trans, int sa=1, int sd=0, dbl s=0.)
	: CBox<Dim>(Trans), StrainAxis(sa), StrainDirection(sd), Strain(s)
{}

template <int Dim>
CLeesEdwardsBox<Dim>::CLeesEdwardsBox(const CLeesEdwardsBox &box) 
	: CBox<Dim>(box), StrainAxis(box.StrainAxis), StrainDirection(box.StrainDirection), Strain(box.Strain)
{}

template <int Dim>	
const CLeesEdwardsBox<Dim>& CLeesEdwardsBox<Dim>::operator=(const CLeesEdwardsBox<Dim> &box)
{
	if(this != &box)
	{
		CBox<Dim>::operator=(box);
		StrainAxis = box.StrainAxis;
		StrainDirection = box.StrainDirection;
		Strain = box.Strain;
	}
	return *this;
}

//functions to read and write box configurations
template <int Dim>
string CLeesEdwardsBox<Dim>::GetName()
{
	stringstream name;
	name << "LeesEdwardsBox_" << Dim << "D";
	return name.str();
}
	
template <int Dim>
string CLeesEdwardsBox<Dim>::DataToString() const
{
	stringstream ss;
	ss << GetName() << ":" << CBox<Dim>::BoxSymmetry; //FIX!!!
	return ss.str();
}
	
template <int Dim>
void CLeesEdwardsBox<Dim>::StringToData(string Data)
{
	vector<string> split = SplitString(Data,":");
	CBox<Dim>::BoxSymmetry = atoi(split[1].c_str()); //FIX!!!
}

template <int Dim>
CLeesEdwardsBox<Dim> *CLeesEdwardsBox<Dim>::Clone() const
{
	return new CLeesEdwardsBox<Dim>( *this );
}


//	class PeriodicBCs
//		This class implements a function to move particle positions
//			according to periodic boundary conditions
//		It is a templated class with a partial specialization
//			to handle the common case of full periodic BCs.
template<int Dim> class PeriodicBCs
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
		assert(Points.cols()%Dim==0);
		int np = Points.cols()/Dim;
		for(int i=0; i<np; i++)
			for(int dd=0; dd<Dim; dd++)
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
template <int Dim>
inline void CLeesEdwardsBox<Dim>::MoveParticles(Eigen::VectorXd &Points,const Eigen::VectorXd &Displacements) const
{
	Points += Displacements;
	for(int i = 0; i< Points.cols() ; i++)
		Points(i) -= floor(Points(i));
}

template <int Dim>
inline void CLeesEdwardsBox<Dim>::MoveParticle(dvecBlock Point, dvec const &Displacement) const
{
	Point += Displacement;
	for(int dd=0; dd<Dim; dd++)
		Point(dd) -= floor(Point(dd));
}




//These two methods assume:
//	1. the box has dimensions 1.
//	2. all particles are within the box.
//
template <int Dim>
inline void CLeesEdwardsBox<Dim>::MinimumDisplacement(const dvec &PointA,const dvec &PointB, dvec &Displacement) const
{
	Displacement = PointA-PointB;
	for(int i = 0; i < Dim ; i++)
	{
		if(i==ShearDirection && abs(Displacement(ShearAxis))>0.5)
		{
			if(PointA(ShearAxis) > 0.5)
				Displacement(i) -= Shear;
			else
				Displacement(i) += Shear;
		}
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
	}
}

template <int Dim>
inline void CLeesEdwardsBox<Dim>::MinimumDisplacement(const dvecBlock &PointA,const dvecBlock &PointB, dvec &Displacement) const
{
	Displacement = PointA-PointB;
	for(int i = 0; i < Dim ; i++)
	{
		if(i==ShearDirection && abs(Displacement(ShearAxis))>0.5)
		{
			if(PointA(ShearAxis) > 0.5)
				Displacement(i) -= Shear;
			else
				Displacement(i) += Shear;
		}
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
	}
//	Displacement = PointA-PointB;
//	for(int i = 0; i < Dim ; i++)
//		if(abs(Displacement(i))>0.5)
//			Displacement(i)-=sgn(Displacement(i));
}

template <int Dim>
void CLeesEdwardsBox<Dim>::GetPeriodicDimensions(std::vector<int> &pdims) const
{
	pdims.clear();
	for(int i=0; i<Dim; ++i)
		pdims.push_back(i);
}

}

#endif //LEES_EDWARDS_BOX
