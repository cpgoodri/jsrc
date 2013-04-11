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

#include "../Resources/std_include.h"
#include "../Resources/MersenneTwister.h"
#include "Potential.h"
#include <map>
#include <string>
#include <sstream>
#include <stdlib.h>

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

class CHarmonicPotential : public CPotential
{	
private:
//the interaction strength
	dbl epsilon;

public:
//constructors and copy operators
	CHarmonicPotential();
	CHarmonicPotential(dbl _e);
	CHarmonicPotential(const CHarmonicPotential &box);
	
	const CHarmonicPotential &operator=(const CHarmonicPotential &box);

//functions to write potential configurations
	string DataToString();
	
//functions to read potential configurations
 	void StringToData(string Data);
 	CPotential *Create();
	
//functions to compute various derivatives
	dbl Compute(dbl dr,dbl r) const;
	dbl ComputeFirstDerivative(dbl dr,dbl r) const;
	dbl ComputeSecondDerivative(dbl dr,dbl r) const;
	dbl ComputeThirdDerivative(dbl dr,dbl r) const;

//functions to compute multiple derivatives at a time
//WARNING: the 4th argument returns f = -d2U/dr2.
	void ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &f) const;
	void ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &f, dbl &k) const;
	
//functions to get other properties of the potential
	dbl ComputeSupport(dbl r) const;

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

//functions to write potential configurations
string CHarmonicPotential::DataToString()
{
	stringstream ss;
	ss << "HarmonicPotential:" << epsilon;
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

void CHarmonicPotential::ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &f) const
{
	dbl delta = 1-dr/r;
	E = 0.5*epsilon*POW2(delta);
	f = epsilon*delta/r;
}

void CHarmonicPotential::ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &f, dbl &k) const
{
	dbl delta = 1-dr/r;
	E = 0.5*epsilon*POW2(delta);
	f = epsilon*delta/r;
	k = epsilon/POW2(r);
}

dbl CHarmonicPotential::ComputeSupport(dbl r) const
{
	return r;
}


}

#endif
