#ifndef STDDATA_DB_H
#define STDDATA_DB_H

#ifndef DONT_USE_NETCDF

#include "Database.h"
#include "../Computers/Data.h"
#include "../Computers/cijkl.h"


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
class CStdDataDatabase : public CDatabase
{
private:
	typedef CStdData<Dim> DATA;
	NcDim *recDim, *dm2Dim, *NcijklDim;
	NcVar *NPpVar, *NcVar, *VolumeVar, *EnergyVar, *PressureVar, *StressVar, *FabricVar, *MaxGradVar, *CijklVar;

public:
	CStdDataDatabase(string fn="temp.nc", NcFile::FileMode mode=NcFile::ReadOnly);

private:
	virtual void SetDimVar();
	virtual void GetDimVar();

public:
	virtual int GetNumRecs() const;

	virtual void Write(DATA const &data, int rec=-1);
	
	virtual void ReadFirst(DATA &data);
	virtual void ReadLast(DATA &data);
	virtual void Read(DATA &data, int rec);
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CStdDataDatabase<Dim>::CStdDataDatabase(string fn, NcFile::FileMode mode)
	: CDatabase(fn,mode)
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
void CStdDataDatabase<Dim>::SetDimVar()
{
	assert(Mode==NcFile::Replace||Mode==NcFile::New);

	//Set the dimensions
	recDim		= File.add_dim("rec");
	dm2Dim		= File.add_dim("dim2", Dim*Dim);
	NcijklDim	= File.add_dim("Ncijkl", cCIJKL<Dim>::num_constants);


	//Set the variables
	NPpVar		= File.add_var("NPp",		ncInt,		recDim);
	NcVar		= File.add_var("Nc",		ncInt,		recDim);
	VolumeVar	= File.add_var("Volume",	ncDouble,	recDim);
	EnergyVar	= File.add_var("Energy",	ncDouble,	recDim);
	PressureVar	= File.add_var("Pressure",	ncDouble,	recDim);
	StressVar	= File.add_var("Stress",	ncDouble,	recDim, dm2Dim);
	FabricVar	= File.add_var("Fabric",	ncDouble,	recDim, dm2Dim);
	MaxGradVar	= File.add_var("MaxGrad",	ncDouble,	recDim);
	CijklVar	= File.add_var("Cijkl",		ncDouble,	recDim, NcijklDim);
}

template <int Dim>
void CStdDataDatabase<Dim>::GetDimVar()
{
	assert(Mode==NcFile::ReadOnly||Mode==NcFile::Write);

	//Get the dimensions
	if(		!(recDim		= File.get_dim("rec"))
		||	!(dm2Dim		= File.get_dim("dim2"))
		||	!(NcijklDim	= File.get_dim("Ncijkl"))
	  ){
		printf("WARNING: Trouble getting the dimensions for the netCDF file.\n");
		FAILFLAG = true;
	}

	//Get the variables
	if(		!(NPpVar		= File.get_var("NPp"))
		||	!(NcVar		= File.get_var("Nc"))
		||	!(VolumeVar	= File.get_var("Volume"))
		||	!(EnergyVar	= File.get_var("Energy"))
		||	!(PressureVar	= File.get_var("Pressure"))
		||	!(StressVar	= File.get_var("Stress"))
		||	!(FabricVar	= File.get_var("Fabric"))
		||	!(MaxGradVar	= File.get_var("MaxGrad"))
		||	!(CijklVar	= File.get_var("Cijkl"))
	  ){
		printf("WARNING: Trouble getting the variables for the netCDF file.\n");
		FAILFLAG = true;
	}
}

template <int Dim>
int CStdDataDatabase<Dim>::GetNumRecs() const
{
	return (int)recDim->size();
}

template <int Dim>
void CStdDataDatabase<Dim>::Write(DATA const &data, int rec)
{
	assert(Mode==NcFile::Replace||Mode==NcFile::Write||Mode==NcFile::New);
	if(rec<0)	rec = recDim->size();

	NPpVar		->put_rec(&data.NPp,			rec);
	NcVar		->put_rec(&data.Nc,				rec);
	VolumeVar	->put_rec(&data.Volume,			rec);
	EnergyVar	->put_rec(&data.Energy,			rec);
	PressureVar	->put_rec(&data.Pressure,		rec);
	StressVar	->put_rec(data.Stress.data(),	rec);
	FabricVar	->put_rec(data.Fabric.data(),	rec);
	MaxGradVar	->put_rec(&data.MaxGrad,		rec);
	CijklVar	->put_rec(&data.cijkl.a[0],		rec);

	File.sync();
}

template <int Dim>
void CStdDataDatabase<Dim>::ReadFirst(DATA &data)
{
	int rec = 0;
	Read(data, rec);
}

template <int Dim>
void CStdDataDatabase<Dim>::ReadLast(DATA &data)
{
	int rec = recDim->size()-1;
	Read(data, rec);
}

template <int Dim>
void CStdDataDatabase<Dim>::Read(DATA &data, int rec)
{
	assert(FAILFLAG == false);
	assert(Mode==NcFile::ReadOnly);

	NPpVar		-> set_cur(rec);
	NcVar		-> set_cur(rec);
	VolumeVar	-> set_cur(rec);
	EnergyVar	-> set_cur(rec);
	PressureVar	-> set_cur(rec);
	StressVar	-> set_cur(rec);
	FabricVar	-> set_cur(rec);
	MaxGradVar	-> set_cur(rec);
	CijklVar	-> set_cur(rec);

	NPpVar		-> get(&data.NPp,			1);
	NcVar		-> get(&data.Nc,			1);
	VolumeVar	-> get(&data.Volume,		1);
	EnergyVar	-> get(&data.Energy,		1);
	PressureVar	-> get(&data.Pressure,		1);
	StressVar	-> get(data.Stress.data(),	1, dm2Dim->size());
	FabricVar	-> get(data.Fabric.data(),	1, dm2Dim->size());
	MaxGradVar	-> get(&data.MaxGrad,		1);
	CijklVar	-> get(&data.cijkl.a[0],	1, NcijklDim->size());

}


}


#endif //DONT_USE_NETCDF

#endif //STDDATA_DB_H

