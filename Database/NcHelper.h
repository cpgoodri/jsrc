#ifndef NCHELPER_H
#define NCHELPER_H

#include <netcdf>

namespace LiuJamming
{
using namespace netCDF;


class NcHelper
{
public:
	static inline NcVar addVar(NcFile &File, string const &name, NcType const &ncType, NcDim const &d1)
	{
		std::vector<NcDim> dv;
		dv.push_back(d1);
		return File.addVar(name, ncType, dv);
	};
	static inline NcVar addVar(NcFile &File, string const &name, NcType const &ncType, NcDim const &d1, NcDim const &d2)
	{
		std::vector<NcDim> dv;
		dv.push_back(d1);
		dv.push_back(d2);
		return File.addVar(name, ncType, dv);
	};
	static inline NcVar addVar(NcFile &File, string const &name, NcType const &ncType, NcDim const &d1, NcDim const &d2, NcDim const &d3)
	{
		std::vector<NcDim> dv;
		dv.push_back(d1);
		dv.push_back(d2);
		dv.push_back(d3);
		return File.addVar(name, ncType, dv);
	};
	static inline NcVar addVar(NcFile &File, string const &name, NcType const &ncType, NcDim const &d1, NcDim const &d2, NcDim const &d3, NcDim const &d4)
	{
		std::vector<NcDim> dv;
		dv.push_back(d1);
		dv.push_back(d2);
		dv.push_back(d3);
		dv.push_back(d4);
		return File.addVar(name, ncType, dv);
	};

	template<typename T>
	static inline void putRec(NcVar &Var, int rec, T const *data)
	{
		std::vector<NcDim> Dims = Var.getDims();
		assert(Dims[0].isUnlimited());

		//Set the start vector
		std::vector<size_t> start(Dims.size(), 0);
		start[0] = rec;

		//Set the count vector
		std::vector<size_t> count;
		count.push_back(1);
		for(int d=1;d<(int)Dims.size(); ++d)
			count.push_back( Dims[d].getSize() );

		//Write the data
		Var.putVar(start, count, data);
	}

	template<typename T>
	static inline void getRec(NcVar &Var, int rec, T *data)
	{
		std::vector<NcDim> Dims = Var.getDims();
		assert(Dims[0].isUnlimited());

		//Set the start vector
		std::vector<size_t> start(Dims.size(), 0);
		start[0] = rec;

		//Set the count vector
		std::vector<size_t> count;
		count.push_back(1);
		for(int d=1;d<(int)Dims.size(); ++d)
			count.push_back( Dims[d].getSize() );

		//Read the data
		Var.getVar(start, count, data);
	}

};







}


#endif //NCHELPER_H
