#include "BasePotential.h"

#ifndef HARMONIC_POTENTIAL
#define HARMONIC_POTENTIAL

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

//! Class to implement harmonic interactions.

class CHarmonicPotential : public CPotential
{	
private:
	dbl epsilon;	//!<the interaction strength

public:
//! @name Constructors and copy operators
///@{
	CHarmonicPotential();
	CHarmonicPotential(dbl _e);
	CHarmonicPotential(const CHarmonicPotential &box);
	const CHarmonicPotential &operator=(const CHarmonicPotential &box);
///@}

//! @name Functions to read and write potential configurations
///@{
	static string GetName();
	string DataToString() const;
 	void StringToData(string Data);
 	CPotential *Create();
	
///@}
	
//! @name Functions to compute various derivatives
///@{ 
	dbl  Compute(dbl dr,dbl r) const;											//!<Compute the 0th derivative
	dbl  ComputeFirstDerivative(dbl dr,dbl r) const;							//!<Compute the 1st derivative
	dbl  ComputeSecondDerivative(dbl dr, dbl r) const;							//!<Compute the 2nd derivative
	dbl  ComputeThirdDerivative(dbl dr, dbl r) const;							//!<Compute the 3th derivative
	void ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const;				//!<Compute the 0th and 1st derivatives
	void ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const;	//!<Compute the 0th and 1st and 2nd derivatives
	
///@}

//! @name Misc.
///@{
	//!Get the maximum distance between two particles of unit diameter such that they interact
	dbl ComputeSupport() const;

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
CHarmonicPotential::CHarmonicPotential()
{
	epsilon = 1.0;
}

CHarmonicPotential::CHarmonicPotential(dbl _e)
{
	epsilon = _e;
}

CHarmonicPotential::CHarmonicPotential(const CHarmonicPotential &pot) : epsilon(pot.epsilon), CPotential(pot)
{
}
	
const CHarmonicPotential &CHarmonicPotential::operator=(const CHarmonicPotential &pot)
{
	epsilon = pot.epsilon;
	CPotential::operator=(pot);
	return *this;
}

string CHarmonicPotential::GetName()
{
	string s = "HarmonicPotential";
	return s;
}

//functions to write potential configurations
string CHarmonicPotential::DataToString() const
{
	stringstream ss;
	ss << GetName() << ":" << epsilon;
	//ss << "HarmonicPotential:" << epsilon;
	return ss.str();
}
	
//functions to read potential configurations
void CHarmonicPotential::StringToData(string Data)
{
	vector<string> split = SplitString(Data,":");
	epsilon = atof(split[1].c_str());
}
 
CPotential *CHarmonicPotential::Create()
{
	return new CHarmonicPotential();
}
	
//functions to compute various derivatives
dbl CHarmonicPotential::Compute(dbl dr,dbl r) const
{
	return 0.5*epsilon*POW2(1-dr/r);
}

dbl CHarmonicPotential::ComputeFirstDerivative(dbl dr,dbl r) const
{
	return -epsilon/r*(1-dr/r);
}

dbl CHarmonicPotential::ComputeSecondDerivative(dbl dr,dbl r) const
{
	return epsilon/POW2(r);
}

dbl CHarmonicPotential::ComputeThirdDerivative(dbl dr,dbl r) const
{
	return 0.;
}

void CHarmonicPotential::ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const
{
	dbl delta = 1-dr/r;
	E = 0.5*epsilon*POW2(delta);
	g = -epsilon*delta/r;
}

void CHarmonicPotential::ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const
{
	dbl delta = 1-dr/r;
	E = 0.5*epsilon*POW2(delta);
	g = -epsilon*delta/r;
	k = epsilon/POW2(r);
}

dbl CHarmonicPotential::ComputeSupport() const
{
	return 1.;
}


}

#endif
