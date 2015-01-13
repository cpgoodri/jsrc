#ifndef HARMONIC_SPRING_POTENTIAL
#define HARMONIC_SPRING_POTENTIAL

//#include "BaseSpringPotential.h"

/////////////////////////////////////////////////////////////////////////////////
//Harmonic Potential class. 
//
//Description
//		This class describes particles that interact via a central force repulsion.
//		It takes a single parameter, the interaction strength \epsilon. The string
//		denoting this potential is "HarmonicPotential"
//
//Global Variables
//
//Variables
//		Interaction strength \epsilon as a dbl
//		
//Implements
//		Calculating various orders of derivatives of the potential
//		Calculating the support of the potential as a function of radius.
//		Convert from the parameter \epsilon to a string reprentation and vice versa
//
//
/////////////////////////////////////////////////////////////////////////////////

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
//
//
//	This implements the two-sided spring potential
//
//	U(dr, dr0, k0) = [ epsilon*k/(alpha*(alpha-1)) ] * (dr0 - dr)^alpha     if dr0>dr
//	U(dr, dr0, k0) = [ epsilon*k/(alpha*(alpha-1)) ] * (dr - dr0)^alpha     if dr0<dr
//
//	with alpha=2. 
//
//
//

//! Class to implement harmonic interactions.

class CHarmonicSpringPotential //: public CPotential
{	
private:
	dbl epsilon;	//!<the interaction strength

public:
//! @name Constructors and copy operators
///@{
	CHarmonicSpringPotential();
	CHarmonicSpringPotential(dbl _e);
	CHarmonicSpringPotential(const CHarmonicSpringPotential &pot);
	const CHarmonicSpringPotential &operator=(const CHarmonicSpringPotential &pot);
///@}

//! @name Functions to read and write potential configurations
///@{
 	virtual CHarmonicSpringPotential *Clone() const;
	
///@}
	
//! @name Functions to compute various derivatives
///@{ 
	virtual dbl  Compute(dbl dr, dbl dr0, dbl k0) const;													//!<Compute the 0th derivative
	virtual dbl  ComputeFirstDerivative(dbl dr, dbl dr0, dbl k0) const;										//!<Compute the 1st derivative
	virtual dbl  ComputeSecondDerivative(dbl dr, dbl dr0, dbl k0) const;									//!<Compute the 2nd derivative
	virtual void ComputeDerivatives01(dbl dr, dbl dr0, dbl k0, dbl &E, dbl &g) const;						//!<Compute the 0th and 1st derivatives
	virtual void ComputeDerivatives012(dbl dr, dbl dr0, dbl k0, dbl &E, dbl &g, dbl &k) const;				//!<Compute the 0th and 1st and 2nd derivatives
	
///@}

//! @name Misc.
///@{

///@}
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//constructors and copy operators
CHarmonicSpringPotential::CHarmonicSpringPotential()
	: epsilon(1.0)
{}

CHarmonicSpringPotential::CHarmonicSpringPotential(dbl _e)
	: epsilon(_e)
{}

CHarmonicSpringPotential::CHarmonicSpringPotential(const CHarmonicSpringPotential &pot) 
	: epsilon(pot.epsilon)
{}
	
const CHarmonicSpringPotential& CHarmonicSpringPotential::operator=(const CHarmonicSpringPotential &pot)
{
	if(this != &pot)
	{
		epsilon = pot.epsilon;
	}
	return *this;
}

/*
string CHarmonicSpringPotential::GetName()
{
	string s = "HarmonicSpringPotential";
	return s;
}
*/

//Clone the potential opbject and return a pointer to the new copy.
CHarmonicSpringPotential* CHarmonicSpringPotential::Clone() const
{
	return new CHarmonicSpringPotential( *this );
}
	
//functions to compute various derivatives
dbl CHarmonicSpringPotential::Compute(dbl dr, dbl dr0, dbl k0) const
{
	return 0.5*epsilon*k0*POW2(dr0-dr);
}

dbl  CHarmonicSpringPotential::ComputeFirstDerivative(dbl dr, dbl dr0, dbl k0) const 
{
	return -epsilon*k0*(dr0-dr);
}

dbl  CHarmonicSpringPotential::ComputeSecondDerivative(dbl dr, dbl dr0, dbl k0) const 
{
	return epsilon*k0;
}

void CHarmonicSpringPotential::ComputeDerivatives01(dbl dr, dbl dr0, dbl k0, dbl &E, dbl &g) const 
{
	E = 0.5*epsilon*k0*POW2(dr0-dr);
	g = -epsilon*k0*(dr0-dr);
}

void CHarmonicSpringPotential::ComputeDerivatives012(dbl dr, dbl dr0, dbl k0, dbl &E, dbl &g, dbl &k) const 
{
	E = 0.5*epsilon*k0*POW2(dr0-dr);
	g = -epsilon*k0*(dr0-dr);
	k = epsilon*k0;
}




}

#endif //HARMONIC_SPRING_POTENTIAL
