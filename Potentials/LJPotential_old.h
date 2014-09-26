#include "BasePotential.h"

#ifndef LJ_REG_POTENTIAL
#define LJ_REG_POTENTIAL

/////////////////////////////////////////////////////////////////////////////////
//Lennard-Jones Potential class. 
//
//Description
//		This class describes particles that interact via a particular form
//		of a Lennard-Jones interaction, standard in the literature. (Kob-Andersen)
//		It takes a single parameter, the interaction strength \epsilon. The string
//		denoting this potential is "LJPotential"
//		
//		It has been specialized to the standard glass-forming set-up (sigma_AB = 0.8 sigma_AA
//		Sigma_BB = 0.88 sigma_AA, epsilon_AB = 1.5 epsilon, epsilon_BB = 0.5 epsilon
//		Cut off at r_ij=2.5 sigma_ij, and shifted to satisfy v(2.5 sigma_ij) = V'(2.5 sigma_ij) V''(2.5 sigma_ij) = 0.
//
//		NOTE: Current implementation specifically requires the ratio of particle sizes to be set like state.SetRadiiBi(f,1.136364) where f is
//			the fraction of small particles...i.e. the precision on the ratio of the sizes should be 1.0/1.136364
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

class CRegLJPotential : public CPotential
{	
private:
	dbl epsilon;	//!<the interaction strength
	dbl rc;			//!<the cutoff distance
	dbl a,b;		//!<constants for offseting the potential

public:
//! @name Constructors and copy operators
///@{
	CRegLJPotential(dbl _rc=2.5, dbl _e=1.0);
	CRegLJPotential(const CRegLJPotential &box);
	const CRegLJPotential &operator=(const CRegLJPotential &box);
///@}

//! @name Functions to read and write potential configurations
///@{
	static string GetName();
	virtual string DataToString() const;
 	virtual void StringToData(string Data);
 	virtual CRegLJPotential *Clone() const;
	
///@}
	
//! @name Functions to compute various derivatives
///@{ 
	dbl  Compute(dbl dr,dbl r) const;											//!<Compute the 0th derivative
	dbl  ComputeFirstDerivative(dbl dr,dbl r) const;							//!<Compute the 1st derivative
	dbl  ComputeSecondDerivative(dbl dr, dbl r) const;							//!<Compute the 2nd derivative
	dbl  ComputeThirdDerivative(dbl dr, dbl r) const;							//!<Compute the 3th derivative
	void ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const;				//!<Compute the 0th and 1st derivatives
	void ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const;	//!<Compute the 0th and 1st and 2nd derivatives
	void ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const;	//!<Compute the 0th and 1st and 2nd and 3rd derivatives

///@}

//! @name Misc.
///@{
	//!Get the maximum distance between two particles of unit diameter such that they interact
	dbl ComputeSupport() const;
	
	//!Set sigma and determine if rad1, rad2 and rlen2 are such that there is overlap
	bool Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const;
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
CRegLJPotential::CRegLJPotential(dbl _rc, dbl _e)
	: epsilon(_e), rc(_rc)
{
	init_ab();
}

CRegLJPotential::CRegLJPotential(const CRegLJPotential &pot) 
	: epsilon(pot.epsilon), 
	  rc(pot.rc),
	  CPotential(pot) //CPG comment: is is necessary to call the CPotential constructor? If so, then this should be done in CHarmonicPotential as well.
{
	init_ab();
}
	
const CRegLJPotential &CRegLJPotential::operator=(const CRegLJPotential &pot)
{
	if(this != &pot)
	{
		rc = pot.rc;
		epsilon = pot.epsilon;
		a = pot.a;
		b = pot.b;
		CPotential::operator=(pot);
	}
	return *this;
}

void CRegLJPotential::init_ab()
{
	dbl rc2 = POW2(rc);
	dbl rc6 = POW3(rc2);
	dbl rc12 = POW2(rc6);

	b = 48./(rc12*rc) - 12./(rc6*rc);
	a = -4.*(1./rc12 - 1./rc6) - b*rc;
}

string CRegLJPotential::GetName()
{
	string s = "LJPotential";
	return s;
}

//functions to write potential configurations
string CRegLJPotential::DataToString() const
{
	stringstream ss; 
	ss << GetName() << STRING_DELIMITER << ConvertDblToHexString(epsilon) << STRING_DELIMITER << ConvertDblToHexString(rc);
	return FillString(ss);
}
	
//functions to read potential configurations
void CRegLJPotential::StringToData(string Data)
{
	vector<string> split = SplitString(Data,":");
	//epsilon = atof(split[1].c_str());
	//CPG comment: I think this might be the problem. Try switching to the following:
	epsilon = ConvertHexStringToDbl(split[1]);
	rc = ConvertHexStringToDbl(split[2]);
}

//Clone the potential opbject and return a pointer to the new copy.
CRegLJPotential* CRegLJPotential::Clone() const
{
        return new CRegLJPotential( *this );
} 
	
//functions to compute various derivatives
dbl CRegLJPotential::Compute(dbl dr,dbl r) const
{
	dbl dri = r/dr;
	dbl dr2i = POW2(dri);
	dbl dr6i = POW3(dr2i);
	dbl dr12i = POW2(dr6i);

	return epsilon*(4.*(dr12i-dr6i) + a + b*dri);
}

dbl CRegLJPotential::ComputeFirstDerivative(dbl dr,dbl r) const
{
	dbl dr2i = POW2(r/dr);
	dbl dr6i = POW3(dr2i);
	dbl dr12i = POW2(dr6i);

	return epsilon*(24.*(dr6i-2.*dr12i)/dr + b/r);
}

dbl CRegLJPotential::ComputeSecondDerivative(dbl dr,dbl r) const
{	
	dbl dr2i = POW2(r/dr);
	dbl dr6i = POW3(dr2i);
	dbl dr12i = POW2(dr6i);

	return epsilon*(624.*dr12i - 168.*dr6i)/POW2(dr);
}

dbl CRegLJPotential::ComputeThirdDerivative(dbl dr,dbl r) const
{
	assert(false);
}

void CRegLJPotential::ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const
{
	dbl dri = r/dr;
	dbl dr2i = POW2(dri);
	dbl dr6i = POW3(dr2i);
	dbl dr12i = POW2(dr6i);

	E = epsilon*(4.*(dr12i-dr6i) + a + b*dri);
	g = epsilon*(24.*(dr6i-2.*dr12i)/dr + b/r);
}

void CRegLJPotential::ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const
{
	dbl dri = r/dr;
	dbl dr2i = POW2(dri);
	dbl dr6i = POW3(dr2i);
	dbl dr12i = POW2(dr6i);

	E = epsilon*(4.*(dr12i-dr6i) + a + b*dri);
	g = epsilon*(24.*(dr6i-2.*dr12i)/dr + b/r);
	k = epsilon*(624.*dr12i - 168.*dr6i)/POW2(dr);
}

void CRegLJPotential::ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const
{
	dbl dri = r/dr;
	dbl dr2i = POW2(dri);
	dbl dr6i = POW3(dr2i);
	dbl dr12i = POW2(dr6i);

	assert(false);
}


dbl CRegLJPotential::ComputeSupport() const
{
	return rc;
}

bool CRegLJPotential::Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const
{
	sigma = rad1 + rad2;
	return (rlen2<sigma*sigma);
}

}

#endif
