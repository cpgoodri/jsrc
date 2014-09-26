#ifndef SS_ARRAY_DB_H
#define SS_ARRAY_DB_H

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

class CSoftSpotArrayDatabase : public CDatabase
{
private:
	int NP;
	NcDim *recDim, *NPDim;
	NcVar *siiVar, *ssArrayVar;

public:
	CSoftSpotArrayDatabase(int np, string fn="temp.nc", NcFile::FileMode mode=NcFile::ReadOnly);

private:
	virtual void SetDimVar();
	virtual void GetDimVar();

public:
	virtual int GetNumRecs() const;

	virtual void Write(int sii, vector<int> const &ssArray, int rec=-1);

	virtual void ReadFirst(int &sii, vector<int> &ssArray);
	virtual void ReadLast(int &sii, vector<int> &ssArray);
	virtual void Read(int &sii, vector<int> &ssArray, int rec);
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

CSoftSpotArrayDatabase::CSoftSpotArrayDatabase(int np, string fn, NcFile::FileMode mode)
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

void CSoftSpotArrayDatabase::SetDimVar()
{
	assert(Mode==NcFile::Replace||Mode==NcFile::New);

	//Set the dimensions
	recDim = File.add_dim("rec");
	NPDim  = File.add_dim("NP",  NP);

	//Set the variables
	siiVar			= File.add_var("sii",		ncInt,		recDim);
	ssArrayVar		= File.add_var("ssArray",	ncShort,	recDim, NPDim);
}

void CSoftSpotArrayDatabase::GetDimVar()
{
	assert(Mode==NcFile::ReadOnly||Mode==NcFile::Write);

	//Get the dimensions
	recDim = File.get_dim("rec");
	NPDim  = File.get_dim("NP");

	//Get the variables
	siiVar			= File.get_var("sii");
	ssArrayVar		= File.get_var("ssArray");
}

int CSoftSpotArrayDatabase::GetNumRecs() const
{
	return (int)recDim->size();
}

void CSoftSpotArrayDatabase::Write(int sii, vector<int> const &ssArray, int rec)
{
	assert(Mode==NcFile::Replace||Mode==NcFile::Write||Mode==NcFile::New);
	assert((int)ssArray.size() == NP);
	if(rec<0)	rec = recDim->size();

	//Create some temporary storage
	vector<short> ssArray_short(NP);
	vector<short>::iterator			s_it = ssArray_short.begin();
	vector<int>  ::const_iterator	i_it = ssArray.begin();
	for( ; s_it != ssArray_short.end(); ++s_it, ++i_it)
		(*s_it) = (short)(*i_it);
	
	//Write all the data
	siiVar		->put_rec(&sii,					rec);
	ssArrayVar	->put_rec(&ssArray_short[0],	rec);

	File.sync();
}


void CSoftSpotArrayDatabase::ReadFirst(int &sii, vector<int> &ssArray)
{
	int rec = 0;
	Read(sii, ssArray, rec);
}

void CSoftSpotArrayDatabase::ReadLast(int &sii, vector<int> &ssArray)
{
	int rec = recDim->size()-1;
	Read(sii, ssArray, rec);
}

void CSoftSpotArrayDatabase::Read(int &sii, vector<int> &ssArray, int rec)
{
	assert(Mode==NcFile::ReadOnly);
	

	//Read the data from the database
	siiVar			-> set_cur(rec);
	ssArrayVar		-> set_cur(rec);

	vector<short> ssArray_short(NPDim->size());

	siiVar			->get(&sii,					1);
	ssArrayVar		->get(&ssArray_short[0],	1, NPDim->size());

	ssArray.clear();
	ssArray.resize(NPDim->size());
	assert(ssArray_short.size() == ssArray.size());

	vector<short>::const_iterator s_it = ssArray_short.begin();
	vector<int>  ::iterator       i_it = ssArray.begin();
	for( ; s_it != ssArray_short.end(); ++s_it, ++i_it)
		(*i_it) = (int)(*s_it);
}




}


#endif //DONT_USE_NETCDF

#endif //SS_ARRAY_DB_H

