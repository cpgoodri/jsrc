#ifndef STATIC_DB_H
#define STATIC_DB_H

#ifndef DONT_USE_NETCDF

#include "Database.h"
#include "../State/StaticState.h"


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
class CStaticDatabase : public CDatabase
{
//private:
public:
	typedef CStaticState<Dim> STATE;
	int NP;
	NcDim *recDim, *dimDim, *dm2Dim, *NPDim, *dofDim, *strDim;
	NcVar *posVar, *radVar, *BoxMatrixVar, *BoxStringVar;//, *PotStringVar;

public:
	CStaticDatabase(int np, string fn="temp.nc", NcFile::FileMode mode=NcFile::ReadOnly);

private:
	virtual void SetDimVar();
	virtual void GetDimVar();

public:
	virtual int GetNumRecs() const;

	virtual void Write(STATE const &s, int rec=-1);

	virtual void ReadFirst(STATE &s);
	virtual void ReadLast(STATE &s);
	virtual void Read(STATE &s, int rec);

	//Read parts of a state
	virtual void ReadBoxMatrix(Eigen::Matrix<dbl,Dim,Dim> &trans, int rec);
	

	static int ExtractSystemSize(std::string fn="temp.nc");

};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CStaticDatabase<Dim>::CStaticDatabase(int np, string fn, NcFile::FileMode mode)
	: CDatabase(fn,mode),
	  NP(np)
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
void CStaticDatabase<Dim>::SetDimVar()
{
	assert(Mode==NcFile::Replace||Mode==NcFile::New);

	//Set the dimensions
	recDim = File.add_dim("rec");
	dimDim = File.add_dim("dim", Dim);
	dm2Dim = File.add_dim("dim2", Dim*Dim);
	NPDim  = File.add_dim("NP",  NP);
	dofDim = File.add_dim("dof", NP*Dim);
	strDim = File.add_dim("StringSize", DB_STRING_SIZE);

	//Set the variables
	posVar			= File.add_var("pos",		ncDouble,	recDim, dofDim);
	radVar			= File.add_var("rad",		ncDouble,	recDim, NPDim );
	BoxMatrixVar	= File.add_var("BoxMatrix",	ncDouble,	recDim, dm2Dim);
	BoxStringVar	= File.add_var("BoxString",	ncChar,		recDim, strDim);
//	PotStringVar	= File.add_var("PotString",	ncChar,		recDim, strDim);
}

template <int Dim>
void CStaticDatabase<Dim>::GetDimVar()
{
	assert(Mode==NcFile::ReadOnly||Mode==NcFile::Write);

	//Get the dimensions
	if(		!(recDim = File.get_dim("rec"))
		||	!(dimDim = File.get_dim("dim"))
		||	!(dm2Dim = File.get_dim("dim2"))
		||	!(NPDim  = File.get_dim("NP"))
		||	!(dofDim = File.get_dim("dof"))
		||	!(strDim = File.get_dim("StringSize"))
	  ){
		printf("WARNING: Trouble getting the dimensions for the netCDF file.\n");
		FAILFLAG = true;
	}

	//Get the variables
	if(		!(posVar			= File.get_var("pos"))
		||	!(radVar			= File.get_var("rad"))
		||	!(BoxMatrixVar	= File.get_var("BoxMatrix"))
		||	!(BoxStringVar	= File.get_var("BoxString"))
	  ){
		printf("WARNING: Trouble getting the variables for the netCDF file.\n");
		FAILFLAG = true;
	}
}

/*
template <int Dim>
void CStaticDatabase<Dim>::GetDimVar()
{
	assert(Mode==NcFile::ReadOnly||Mode==NcFile::Write);

	//Get the dimensions
	recDim = File.get_dim("rec");
	dimDim = File.get_dim("dim");
	dm2Dim = File.get_dim("dim2");
	NPDim  = File.get_dim("NP");
	dofDim = File.get_dim("dof");
	strDim = File.get_dim("StringSize");

	//Get the variables
	posVar			= File.get_var("pos");
	radVar			= File.get_var("rad");
	BoxMatrixVar	= File.get_var("BoxMatrix");
	BoxStringVar	= File.get_var("BoxString");
}
*/

template <int Dim>
int CStaticDatabase<Dim>::GetNumRecs() const
{
	return (int)recDim->size();
}

static void CopyString(char *cstr, string &str, int max_size, char background)
{
	int i;
	for(i=0; i<str.size()&&i<max_size; ++i)
		cstr[i] = str[i];
	for( ; i<max_size; ++i)
		cstr[i] = background;
}

template <int Dim>
void CStaticDatabase<Dim>::Write(STATE const &s, int rec)
{
	assert(Mode==NcFile::Replace||Mode==NcFile::Write||Mode==NcFile::New);
	assert(s.GetParticleNumber() == NP);
	if(rec<0)	rec = recDim->size();

	//Create some temporary storage
	Eigen::Matrix<dbl,Dim,Dim> trans;
	char BoxCString[DB_STRING_SIZE];
	string BoxString;

	//Prepare the box and potential data
	//BoxMatrix
	s.GetBox()->GetTransformation(trans);
	//Box String
	BoxString = s.GetBox()->DataToString();
	CopyString(BoxCString, BoxString, DB_STRING_SIZE, ':');

	//Write all the data
	posVar		->put_rec(&s.Positions[0],	rec);
	radVar		->put_rec(&s.Radii[0],		rec);
	BoxMatrixVar->put_rec(trans.data(),		rec);
	BoxStringVar->put_rec(&BoxCString[0],	rec);
	s.GetPotential()->NetCDFWrite(File, rec);

	File.sync();
}


template <int Dim>
void CStaticDatabase<Dim>::ReadFirst(STATE &s)
{
	int rec = 0;
	Read(s, rec);
}

template <int Dim>
void CStaticDatabase<Dim>::ReadLast(STATE &s)
{
	int rec = recDim->size()-1;
	Read(s, rec);
}

template <int Dim>
void CStaticDatabase<Dim>::Read(STATE &s, int rec)
{
	assert(FAILFLAG == false);
	assert(Mode==NcFile::ReadOnly);
	assert(s.GetParticleNumber() == NP);
	
	//Check Dimensions
	//	todo
	
	Eigen::Matrix<dbl,Dim,Dim> trans;
	char BoxCString[DB_STRING_SIZE];

	//Read the data from the database
	posVar			-> set_cur(rec);
	radVar			-> set_cur(rec);
	BoxMatrixVar	-> set_cur(rec);
	BoxStringVar	-> set_cur(rec);

	posVar			->get(&s.Positions[0],	1, dofDim->size());
	radVar			->get(&s.Radii[0],		1, NPDim->size());
	BoxMatrixVar	->get(trans.data(),		1, dm2Dim->size());
	BoxStringVar	->get(&BoxCString[0],	1, strDim->size());

	//Set the box
	s.Box = CBox<Dim>::SetFromStringAndMatrix(BoxCString,trans);
	s.Box->SetSymmetry(CBox<Dim>::QUADRILATERAL); //This sometimes isn't taken care of by the constructors... not sure why...
//	printf("111: %i\n", s.Box->GetSymmetry());

	//Set the potential
	s.Potential = CPotential::NetCDFRead(File, rec);
}

template <int Dim>
void CStaticDatabase<Dim>::ReadBoxMatrix(Eigen::Matrix<dbl,Dim,Dim> &trans, int rec)
{
	assert(Mode==NcFile::ReadOnly);
	BoxMatrixVar	-> set_cur(rec);
	BoxMatrixVar	->get(trans.data(),		1, dm2Dim->size());
}

	
template <int Dim>
int CStaticDatabase<Dim>::ExtractSystemSize(std::string fn)
{
	NcFile File(fn.c_str(), NcFile::ReadOnly);
	NcDim *NPDim = File.get_dim("NP");
	int NP_temp = (int)NPDim->size();
	return NP_temp;
}










}


#endif //DONT_USE_NETCDF

#endif //STATIC_DB_H

