#ifndef STATIC_DB_H
#define STATIC_DB_H

#ifndef DONT_USE_NETCDF

#include "Database.h"


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
private:
	typedef CStaticState<Dim> STATE;
	int NP;
	NcDim recDim, dimDim, dm2Dim, NPDim, dofDim, strDim;
	NcVar posVar, radVar, BoxMatrixVar, BoxStringVar, PotStringVar;

	int Current;


public:
	CStaticDatabase(int np, string fn="temp.nc", NcFile::FileMode mode=NcFile::read, NcFile::FileFormat format=NcFile::nc4);

private:
	void SetDimVar();
	void GetDimVar();

public:
	void SetCurrentRec(int r);
	int  GetCurrentRec();

	void WriteState(STATE const &c, int rec=-1);
	void ReadState(STATE &c, int rec);
	void ReadNextState(STATE &c);
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CStaticDatabase<Dim>::CStaticDatabase(int np, string fn, NcFile::FileMode mode, NcFile::FileFormat format)
	: CDatabase(fn,mode,format),
	  NP(np),
	  Current(0)
{
	switch(Mode)
	{
		case _read:
		case _write:
			GetDimVar();
			break;
		case _replace:
		case _newFile:
			SetDimVar();
			break;
		default:
			assert(false);
	}
};

template <int Dim>
void CStaticDatabase<Dim>::SetDimVar()
{
	try{
		assert(Mode==_replace||Mode==_newFile);

		//Set the dimensions
		recDim = File.addDim("rec");
		dimDim = File.addDim("dim", Dim);
		dm2Dim = File.addDim("dim2", Dim*Dim);
		NPDim  = File.addDim("NP",  NP);
		dofDim = File.addDim("dof", NP*Dim);
		strDim = File.addDim("StringSize", DB_STRING_SIZE);

		//Set the variables
		posVar			= NcHelper::addVar(File, "pos",			ncDouble,	recDim, dofDim);
		radVar			= NcHelper::addVar(File, "rad",			ncDouble,	recDim, NPDim );
		BoxMatrixVar	= NcHelper::addVar(File, "BoxMatrix",	ncDouble,	recDim, dm2Dim);
		BoxStringVar	= NcHelper::addVar(File, "BoxString",	ncChar,		recDim, strDim);
		PotStringVar	= NcHelper::addVar(File, "PotString",	ncChar,		recDim, strDim);
	}catch (NcException& e){    e.what();   }
}

template <int Dim>
void CStaticDatabase<Dim>::GetDimVar()
{
	try{
		assert(Mode==_read||Mode==_write);

		//Get the dimensions
		recDim = File.getDim("rec");
		dimDim = File.getDim("dim");
		dm2Dim = File.getDim("dim2");
		NPDim  = File.getDim("NP");
		dofDim = File.getDim("dof");
		strDim = File.getDim("StringSize");

		//Get the variables
		posVar			= File.getVar("pos");
		radVar			= File.getVar("rad");
		BoxMatrixVar	= File.getVar("BoxMatrix");
		BoxStringVar	= File.getVar("BoxString");
		PotStringVar	= File.getVar("PotString");
	}catch (NcException& e){    e.what();   }
}


template <int Dim>
void CStaticDatabase<Dim>::SetCurrentRec(int r)
{
	Current = r;
}

template <int Dim>
int CStaticDatabase<Dim>::GetCurrentRec()
{
	return Current;
}

void CopyString(char *cstr, string &str, int max_size, char background)
{
	int i;
	for(i=0; i<str.size()&&i<max_size; ++i)
		cstr[i] = str[i];
	for( ; i<max_size; ++i)
		cstr[i] = background;
}

template <int Dim>
void CStaticDatabase<Dim>::WriteState(STATE const &s, int rec)
{
	try{
		assert(Mode==_replace||Mode==_write||Mode==_newFile);
		assert(s.GetParticleNumber() == NP);
		if(rec<0)	rec = recDim.getSize();

		//Create some temporary storage
		Eigen::Matrix<dbl,Dim,Dim> trans;
		char BoxCString[DB_STRING_SIZE];
		char PotCString[DB_STRING_SIZE];
		string BoxString, PotString;

		//Prepare the box and potential data
		//BoxMatrix
		s.GetBox()->GetTransformation(trans);
		//Box String
		BoxString = s.GetBox()->DataToString();
		CopyString(BoxCString, BoxString, DB_STRING_SIZE, ':');
		//Potential String
		PotString = s.GetPotential()->DataToString();
		CopyString(PotCString, PotString, DB_STRING_SIZE, ':');

		//Write all the data
		NcHelper::putRec(posVar,		rec, &s.Positions[0]);
		NcHelper::putRec(radVar,		rec, &s.Radii[0]);
		NcHelper::putRec(BoxMatrixVar,	rec, trans.data());
		NcHelper::putRec(BoxStringVar,	rec, &BoxCString[0]);
		NcHelper::putRec(PotStringVar,	rec, &PotCString[0]);
	}catch (NcException& e){    e.what();   }
}

template <int Dim>
void CStaticDatabase<Dim>::ReadState(STATE &s, int rec)
{
	try{
		assert(Mode==_read);
		assert(s.GetParticleNumber() == NP);
		
		//Check Dimensions
		
		Eigen::Matrix<dbl,Dim,Dim> trans;
		char BoxCString[DB_STRING_SIZE];
		char PotCString[DB_STRING_SIZE];

		//Read the data from the database
		NcHelper::getRec(posVar,		rec, &s.Positions[0]);
		NcHelper::getRec(radVar,		rec, &s.Radii[0]);
		NcHelper::getRec(BoxMatrixVar,	rec, trans.data());
		NcHelper::getRec(BoxStringVar,	rec, &BoxCString[0]);
		NcHelper::getRec(PotStringVar,	rec, &PotCString[0]);

		//Set the box
		s.Box = CBox<Dim>::SetFromStringAndMatrix(BoxCString,trans);

		//Set the potential
		s.Potential = CPotential::SetFromString(PotCString);
	}catch (NcException& e){    e.what();   }
}




}


#endif //DONT_USE_NETCDF

#endif //STATIC_DB_H

