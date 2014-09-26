#include "BasePotential.h"

#ifndef LJ_REG_POTENTIAL
#define LJ_REG_POTENTIAL

/////////////////////////////////////////////////////////////////////////////////
//Lennard-Jones Potential class. 
//
//Description
//		This class describes particles that interact via a Lennard-Jones interaction.
//		It takes two parameters, the interaction strength \epsilon, and the cutoff distance r_c. 
//		The string denoting this potential is "RegLJPotential"
//		
//		This follows the lj/smooth/linear potential used in LAMMPS, except that \sigma
//		is always determined from the sum of the "radii". This is therefore not easily 
//		extended to the Kob-Anderson interaction, but does allow for polydispersity.
//
//		The potential is given by
//
//		E(r, sigma) = phi(r,sigma) - phi(rc,sigma) - (r-rc) dphi/dr(evaluated at r=rc)
//
//		where
//
//		phi(r,sigma) = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
//
//		This cal also be written as 
//
//		E(r,sigma) = epsilon ( 4[ (sigma/r)^12 - (sigma/4)^6 ] + a + b r/sigma
//
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
	dbl epsilon;	//!<the interaction strength. Default = 1.0
	dbl rc;			//!<the cutoff distance, in units of sigma. Default = 2.5
	dbl a,b;		//!<constants for offseting the potential
	void init_ab();	//!<Initialize a and b

public:
//! @name Constructors and copy operators
///@{
	CRegLJPotential(dbl _rc=2.5, dbl _e=1.0);
	CRegLJPotential(const CRegLJPotential &pot);
	const CRegLJPotential &operator=(const CRegLJPotential &pot);
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
	dbl rci   = 1./rc;
	dbl rc2i  = POW2(rci);
	dbl rc6i  = POW3(rc2i);
	dbl rc12i = POW2(rc6i);

	a = -52.*rc12i + 28.*rc6i;
	b = rci*(48.*rc12i - 24.*rc6i);
}

string CRegLJPotential::GetName()
{
	string s = "RegLJPotential";
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

// dr = r
// r  = sigma
//functions to compute various derivatives
dbl CRegLJPotential::Compute(dbl dr,dbl r) const
{
	dbl rs    = dr/r; // in better notation, this is r/sigma
	dbl rs2i  = 1./POW2(rs);
	dbl rs6i  = POW3(rs2i);
	dbl rs12i = POW2(rs6i);

	return epsilon*(4.*(rs12i-rs6i) + a + b*rs);
}

dbl CRegLJPotential::ComputeFirstDerivative(dbl dr,dbl r) const
{
	dbl rs2i  = POW2(r/dr);
	dbl rs6i  = POW3(rs2i);
	dbl rs12i = POW2(rs6i);

	return epsilon*(24.*(rs6i-2.*rs12i)/dr + b/r);
}

dbl CRegLJPotential::ComputeSecondDerivative(dbl dr,dbl r) const
{	
	dbl rs2i  = POW2(r/dr);
	dbl rs6i  = POW3(rs2i);
	dbl rs12i = POW2(rs6i);

	return epsilon*(624.*rs12i - 168.*rs6i)/POW2(dr);
}

dbl CRegLJPotential::ComputeThirdDerivative(dbl dr,dbl r) const
{
	dbl rs2i  = POW2(r/dr);
	dbl rs6i  = POW3(rs2i);
	dbl rs12i = POW2(rs6i);

	return epsilon*(8736.*rs12i - 1344.*rs6i)/POW3(dr);
}

void CRegLJPotential::ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const
{
	dbl rs    = dr/r; // in better notation, this is r/sigma
	dbl rs2i  = 1./POW2(rs);
	dbl rs6i  = POW3(rs2i);
	dbl rs12i = POW2(rs6i);

	E = epsilon*(4.*(rs12i-rs6i) + a + b*rs);
	g = epsilon*(24.*(rs6i-2.*rs12i)/dr + b/r);
	//k = epsilon*(624.*rs12i - 168.*rs6i)/POW2(dr);
	//t = epsilon*(8736.*rs12i - 1344.*rs6i)/POW3(dr);
}

void CRegLJPotential::ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const
{
	dbl rs    = dr/r; // in better notation, this is r/sigma
	dbl rs2i  = 1./POW2(rs);
	dbl rs6i  = POW3(rs2i);
	dbl rs12i = POW2(rs6i);

	E = epsilon*(4.*(rs12i-rs6i) + a + b*rs);
	g = epsilon*(24.*(rs6i-2.*rs12i)/dr + b/r);
	k = epsilon*(624.*rs12i - 168.*rs6i)/POW2(dr);
	//t = epsilon*(8736.*rs12i - 1344.*rs6i)/POW3(dr);
}

void CRegLJPotential::ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const
{
	dbl rs    = dr/r; // in better notation, this is r/sigma
	dbl rs2i  = 1./POW2(rs);
	dbl rs6i  = POW3(rs2i);
	dbl rs12i = POW2(rs6i);

	E = epsilon*(4.*(rs12i-rs6i) + a + b*rs);
	g = epsilon*(24.*(rs6i-2.*rs12i)/dr + b/r);
	k = epsilon*(624.*rs12i - 168.*rs6i)/POW2(dr);
	t = epsilon*(8736.*rs12i - 1344.*rs6i)/POW3(dr);
}


dbl CRegLJPotential::ComputeSupport() const
{
	return rc;
}

bool CRegLJPotential::Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const
{
	sigma = rad1 + rad2;
	return (rlen2<POW2(sigma*rc));
}

}

#endif
