#ifndef NETWORK_DB_H
#define NETWORK_DB_H

#ifndef DONT_USE_NETCDF

#include "Database.h"
#include "../State/NetworkState.h"


namespace LiuJamming
{



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
class CNetworkDatabase : public CDatabase
{
private:
	typedef CNetworkState<Dim> STATE;
	int NP;
	int NBonds;
	NcDim *recDim, *dimDim, *dm2Dim, *NPDim, *dofDim, *NBDim, *strDim;
	NcVar *posVar, *bondiVar, *bondjVar, *stiffnessVar, *elengthVar, *BoxMatrixVar, *BoxStringVar;//, *PotStringVar;

public:
	CNetworkDatabase(int np, int nbonds, string fn="temp.nc", NcFile::FileMode mode=NcFile::ReadOnly);

private:
	virtual void SetDimVar();
	virtual void GetDimVar();

public:
	virtual int GetNumRecs() const;

	virtual void Write(STATE const &ns, int rec=-1);

	virtual void ReadFirst(STATE &ns);
	virtual void ReadLast(STATE &ns);
	virtual void Read(STATE &ns, int rec);

	//Read parts of a state
	virtual void ReadBoxMatrix(Eigen::Matrix<dbl,Dim,Dim> &trans, int rec);
	

//	static int ExtractSystemSize(std::string fn="temp.nc");
	static void ExtractSystemSize(std::string fn, int &N, int &NBonds);

};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CNetworkDatabase<Dim>::CNetworkDatabase(int np, int nbonds, string fn, NcFile::FileMode mode)
	: CDatabase(fn,mode),
	  NP(np),
	  NBonds(nbonds)
{
	switch(Mode)
	{
		case NcFile::ReadOnly:
		case NcFile::Write:
			GetDimVar();
			break;
		case NcFile::Replace:
		case NcFile::New:
			SetDimVar();
			break;
		default:
			assert(false);
	}
};

template <int Dim>
void CNetworkDatabase<Dim>::SetDimVar()
{
	assert(Mode==NcFile::Replace||Mode==NcFile::New);

	//Set the dimensions
	recDim		= File.add_dim("rec");
	dimDim		= File.add_dim("dim", Dim);
	dm2Dim		= File.add_dim("dim2", Dim*Dim);
	NPDim		= File.add_dim("NP",  NP);
	dofDim		= File.add_dim("dof", NP*Dim);
	NBDim		= File.add_dim("NBonds", NBonds);
	strDim		= File.add_dim("StringSize", DB_STRING_SIZE);

	//Set the variables
	posVar			= File.add_var("pos",		ncDouble,	recDim, dofDim);
	bondiVar		= File.add_var("bondi",		ncInt,		recDim, NBDim);
	bondjVar		= File.add_var("bondj",		ncInt,		recDim, NBDim);
	stiffnessVar	= File.add_var("stiffness",	ncDouble,	recDim, NBDim);
	elengthVar		= File.add_var("elength",	ncDouble,	recDim, NBDim);
	BoxMatrixVar	= File.add_var("BoxMatrix",	ncDouble,	recDim, dm2Dim);
	BoxStringVar	= File.add_var("BoxString",	ncChar,		recDim, strDim);
}

template <int Dim>
void CNetworkDatabase<Dim>::GetDimVar()
{
	assert(Mode==NcFile::ReadOnly||Mode==NcFile::Write);

	//Get the dimensions
	recDim = File.get_dim("rec");
	dimDim = File.get_dim("dim");
	dm2Dim = File.get_dim("dim2");
	NPDim  = File.get_dim("NP");
	dofDim = File.get_dim("dof");
	NBDim  = File.get_dim("NBonds");
	strDim = File.get_dim("StringSize");

	//Get the variables
	posVar			= File.get_var("pos");
	bondiVar		= File.get_var("bondi");
	bondjVar		= File.get_var("bondj");
	stiffnessVar	= File.get_var("stiffness");
	elengthVar		= File.get_var("elength");
	BoxMatrixVar	= File.get_var("BoxMatrix");
	BoxStringVar	= File.get_var("BoxString");
}

template <int Dim>
int CNetworkDatabase<Dim>::GetNumRecs() const
{
	return (int)recDim->size();
}

/*
 *	This is already defined in the CStaticDatabase.h file... its not pretty but I'm just going to comment it out here.
 *
static void CopyString(char *cstr, string &str, int max_size, char background)
{
	int i;
	for(i=0; i<str.size()&&i<max_size; ++i)
		cstr[i] = str[i];
	for( ; i<max_size; ++i)
		cstr[i] = background;
}
*/

template <int Dim>
void CNetworkDatabase<Dim>::Write(STATE const &ns, int rec)
{
	assert(Mode==NcFile::Replace||Mode==NcFile::Write||Mode==NcFile::New);
	assert(ns.GetNodeNumber() == NP);
	assert(ns.GetNBonds() == NBonds);
	if(rec<0)	rec = recDim->size();

	//Create some temporary storage
	Eigen::Matrix<dbl,Dim,Dim> trans;
	char BoxCString[DB_STRING_SIZE];
	string BoxString;

	//Prepare the box and potential data
	//BoxMatrix
	ns.GetBox()->GetTransformation(trans);
	//Box String
	BoxString = ns.GetBox()->DataToString();
	CopyString(BoxCString, BoxString, DB_STRING_SIZE, ':');

	//Write all the data
	posVar			->put_rec(&ns.Positions[0],		rec);
	bondiVar		->put_rec(&ns.Bondi[0],			rec);
	bondjVar		->put_rec(&ns.Bondj[0],			rec);
	stiffnessVar	->put_rec(&ns.Stiffnesses[0],	rec);
	elengthVar		->put_rec(&ns.ELengths[0],		rec);

	BoxMatrixVar->put_rec(trans.data(),		rec);
	BoxStringVar->put_rec(&BoxCString[0],	rec);
//	s.GetPotential()->NetCDFWrite(File,		rec);

	File.sync();
}


template <int Dim>
void CNetworkDatabase<Dim>::ReadFirst(STATE &ns)
{
	int rec = 0;
	Read(ns, rec);
}

template <int Dim>
void CNetworkDatabase<Dim>::ReadLast(STATE &ns)
{
	int rec = recDim->size()-1;
	Read(ns, rec);
}

template <int Dim>
void CNetworkDatabase<Dim>::Read(STATE &ns, int rec)
{
//	assert(Mode==NcFile::ReadOnly);
//	assert(ns.GetNodeNumber() == NP);
	
	//Check Dimensions
	//	todo

	Eigen::Matrix<dbl,Dim,Dim> trans;
	char BoxCString[DB_STRING_SIZE];

	//Read the data from the database
	posVar			-> set_cur(rec);
	bondiVar		-> set_cur(rec);
	bondjVar		-> set_cur(rec);
	stiffnessVar	-> set_cur(rec);
	elengthVar		-> set_cur(rec);
	BoxMatrixVar	-> set_cur(rec);
	BoxStringVar	-> set_cur(rec);

	ns.N = NP;
	ns.Nbonds = NBonds;

	ns.Positions.resize(Dim*NP);
	ns.Bondi.resize(NBonds);
	ns.Bondj.resize(NBonds);
	ns.Stiffnesses.resize(NBonds);
	ns.ELengths.resize(NBonds);

	posVar			->get(&ns.Positions[0],		1, dofDim->size());
	bondiVar		->get(&ns.Bondi[0],			1, NBDim->size());
	bondjVar		->get(&ns.Bondj[0],			1, NBDim->size());
	stiffnessVar	->get(&ns.Stiffnesses[0],	1, NBDim->size());
	elengthVar		->get(&ns.ELengths[0],		1, NBDim->size());
	BoxMatrixVar	->get(trans.data(),			1, dm2Dim->size());
	BoxStringVar	->get(&BoxCString[0],		1, strDim->size());

	//Set the box
	ns.Box = CBox<Dim>::SetFromStringAndMatrix(BoxCString,trans);
	ns.Box->SetSymmetry(CBox<Dim>::QUADRILATERAL); //This sometimes isn't taken care of by the constructors... not sure why...

	if(ns.Potential!=NULL) delete ns.Potential;
	ns.Potential   = new CHarmonicSpringPotential(); //WARNING!! NOT GENERAL. This is the only potential I have implemented so far for springs

//	printf("111: %i\n", s.Box->GetSymmetry());

//	//Set the potential
//	s.Potential = CPotential::NetCDFRead(File, rec);
}

template <int Dim>
void CNetworkDatabase<Dim>::ReadBoxMatrix(Eigen::Matrix<dbl,Dim,Dim> &trans, int rec)
{
	assert(Mode==NcFile::ReadOnly);
	BoxMatrixVar	-> set_cur(rec);
	BoxMatrixVar	->get(trans.data(),		1, dm2Dim->size());
}

	
template <int Dim>
void CNetworkDatabase<Dim>::ExtractSystemSize(std::string fn, int &N, int &NBonds)
{
	NcFile File(fn.c_str(), NcFile::ReadOnly);
	NcDim *NPDim = File.get_dim("NP");
	NcDim *NBDim = File.get_dim("NBonds");
	N = (int)NPDim->size();
	NBonds = (int)NBDim->size();
}










}


#endif //DONT_USE_NETCDF

#endif //NETWORK_DB_H

