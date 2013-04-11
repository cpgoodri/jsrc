#ifndef POTENTIAL

#define POTENTIAL

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
#include "../Resources/MersenneTwister.h"
#include "netcdfcpp.h"
#include <map>
#include <string>
#include <vector>

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

class CPotential
{
private:
//global variables to read box configurations
	static std::map<string,CPotential*> PotentialTypes;

//Functions to check that the correct dimensions, variables, and attributes are present
//in a netCDF file.
	static bool CheckNetCDF(const NcFile &file);
	static void PopulateNetCDF(NcFile &file);
	
public:
//global functions to read box configurations
	static CPotential *Read(const NcFile &file,int record);
	static void AddPotentialType(string type,CPotential *pot);

private:
	
public:
//constructors and copy operators

//functions to write potential configurations
	void Write(NcFile &file,int record); 
	virtual string DataToString() = 0;
	
//functions to read potential configurations
 	virtual void StringToData(string Data) = 0;
 	virtual CPotential *Create() = 0;
	
//functions to compute various derivatives
	virtual dbl Compute(dbl dr,dbl r) const = 0;
	virtual dbl ComputeFirstDerivative(dbl dr,dbl r) const = 0;
	virtual dbl ComputeSecondDerivative(dbl dr, dbl r) const = 0;
	virtual dbl ComputeThirdDerivative(dbl dr, dbl r) const = 0;

//functions to compute multiple derivatives at a time
	virtual void ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const = 0;
	virtual void ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const = 0;
	
//functions to get other properties of the potential
	virtual dbl ComputeSupport(dbl r) const = 0;
	
};

std::map<string,CPotential*> CPotential::PotentialTypes;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


bool CPotential::CheckNetCDF(const NcFile &file)
{
	return (file.get_att("Potential_Populated")!=NULL);
}

void CPotential::PopulateNetCDF(NcFile &file)
{
	file.add_att("Potential_Populated",1);
	
	NcDim *recorddim = file.get_dim("Records");
	NcDim *datadim = file.get_dim("Data_Size");
	if(datadim==NULL)
		datadim = file.add_dim("Data_Size",256);
	
	file.add_var("Potential_Data",ncChar,recorddim,datadim);
}
	
//global functions to read box configurations
CPotential *CPotential::Read(const NcFile &file,int record)
{
	if(!CheckNetCDF(file))
		throw(CException("CPotential<Dim>::ReadPotential","Attempting to read a box from a file that has no appropriate box data."));
	
	if(record>=file.get_dim("Records")->size())
		throw(CException("CPotential<Dim>::ReadPotential","Attempting to read a box from a record that does not exist."));


	NcVar * pot_var = file.get_var("Potential_Data");
	pot_var->set_cur(record);
	char t_str[256];
	pot_var->get(t_str,1,256);
	string str = t_str;
	vector<string> split = SplitString(str,":");
	CPotential *pot = PotentialTypes[split[0]]->Create();
	pot->StringToData(str);
	return pot;
}

void CPotential::AddPotentialType(string type,CPotential *box)
{
	CPotential::PotentialTypes[type] = box;
}


//functions to write potential configurations
void CPotential::Write(NcFile &file,int record)
{
	if(!CheckNetCDF(file))
		PopulateNetCDF(file);
	
	if(record>file.get_dim("Records")->size())
		throw(CException("CPotential<Dim>::ReadPotential","Attempting to read a box from a record that does not exist."));
 
 	NcVar * pot_var = file.get_var("Potential_Data");
	pot_var->set_cur(record);
	pot_var->put(DataToString().c_str(),1,DataToString().length());
}


}

#endif
