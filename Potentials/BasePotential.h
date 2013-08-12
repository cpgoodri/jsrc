#ifndef BASE_POTENTIAL_H
#define BASE_POTENTIAL_H

/////////////////////////////////////////////////////////////////////////////////
//Potential class. 
//
//Description
//		This is a virtual class that describes a generic isotropic potential. 
//		This class implements functions to read and write the potential and its
//		parameters formatted as a string to a NetCDF file. This class provides
//		functions to compute various derivated of the potential with respect 
//		to the distance dr between two particles. This class also has a function
//		that gives the support of the potential. It is generally assumed that all
//		potentials will have compact support.
//
//		Generically the energy will be related to the ratio of the distance between
//		particles to the sum of their radii.
//
//Global Variables
//		A map from string to potential type as an STD map.
//
//Variables
//		
//Implements
//		Reading/Writing to netcdf files.
//
//Virtual Functions
//		Calculating various orders of derivatives of the potential
//		Calculating the support of the potential as a function of radius.
//
//File Formats
//		NetCDF
//			## Note: To denote whether the NetCDF file has been populated with variables,
//			## an attribute "Potential_Populated" gets created.
//			-> One dimensions: records (unlimited).
//			-> One attribute Potential_Populated.
//			-> String of parameters as a variable.
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include "../Resources/Exception.h"
#include "../Resources/Resources.h"
//#include "netcdfcpp.h"

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

//!Abstract base class for pair-potentials.
class CPotential
{
private:
	static std::map<string,CPotential*> PotentialTypes; //!<Global map to reference potentials by name

public:
//! @name Global functions to read and write
///@{
	static CPotential* SetFromString(string str);	//!<From an input string, return a pointer to a newly created potential object with the given parameter values
///@}

//! @name Functions for setting, copying and saving
///@{
	virtual string DataToString() const = 0;		//!<Return a string that contains the data necessary for the potential object (e.g. interaction strength, etc.)
 	virtual void StringToData(string Data) = 0;		//!<From an input string, set the data variables.
 	virtual CPotential* Clone() const = 0;			//!<Return a pointer to a newly created clone of the object.

///@}

//! @name Functions to compute various derivatives
///@{
	virtual dbl  Compute(dbl dr,dbl r) const = 0;													//!<Compute the 0th derivative
	virtual dbl  ComputeFirstDerivative(dbl dr,dbl r) const = 0;									//!<Compute the 1st derivative
	virtual dbl  ComputeSecondDerivative(dbl dr, dbl r) const = 0;									//!<Compute the 2nd derivative
	virtual dbl  ComputeThirdDerivative(dbl dr, dbl r) const = 0;									//!<Compute the 3th derivative
	virtual void ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const = 0;						//!<Compute the 0th and 1st derivatives
	virtual void ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const = 0;			//!<Compute the 0th and 1st and 2nd derivatives
	virtual void ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const = 0;	//!<Compute the 0th and 1st and 2nd and 3rd derivatives
	
///@}

//! @name Misc.
///@{
	//!Get the maximum distance between two particles of unit diameter such that they interact
	virtual dbl ComputeSupport() const = 0;
	
	//!Set sigma and determine if rad1, rad2 and rlen2 are such that there is overlap
	virtual bool Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const = 0;
///@}
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


CPotential* CPotential::SetFromString(string str)
{
	vector<string> split = SplitString(str,":");
	CPotential *pot = PotentialTypes[split[0]]->Clone();
	pot->StringToData(str);
	return pot;
}


}

#endif //BASE_POTENTIAL_H
