#ifndef BOX

#define BOX

/////////////////////////////////////////////////////////////////////////////////
//Box class. 
//
//Description
//		This is a virtual class that describes box shape. The shape of the box is
//		specified by a dxd matrix that describes the mapping from the unit box [0,1]x...x[0,1]
//		to a corresponding parallelpiped. Overloaded classes will specify the box topology
//		by specifying how the various edges of the box connect to one another. This class
//		implements functions to take points from the base space to the image space. Virtual
//		functions are provided to move points and compute minimum displacements between 
//		points. 
//
//		A framework to generically read and write box configurations to and from netcdf
//		files is also implemented. New box objects that inherit from the box class need
//		to register a string with their name to the map of box types. They must also specify
//		how to read/write any parameters (which must all be doubles) to a string. Inherited
//		objects must further implement several copy constructors. Finally, inherited classes must 
//		implement a function to create a new object of the given class.
//
//Global Variables
//		A map from string to box type as an STL map.
//
//Variables
// 		Transformation as a dxd matrix.
//		Inverse transformation as a dxd matrix.
//		
//Implements
//		Reading/Writing to netcdf files.
//		Mapping to and from the unit box to a transformed parallelopiped.
//
//Virtual Functions
//		Calculating minimal distances between points.
//		Moving particles while respecting the boundary conditions.
//		Mapping the system parameters to and from a string.
//
//File Formats
//		NetCDF
//			## Note: only files of a single Dimension may be stored in a given 
//			## file. To denote whether the NetCDF file has been populated with variables,
//			## dimensions, etc... an attribute "Box_Populated" gets created.
//			-> Two dimensions: records (unlimited) and dimension.
//			-> One attribute Box_Populated.
//			-> Transformation as a variable
//			-> String of parameters as a string.
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include <Eigen/LU>
#include "../Resources/MersenneTwister.h"
#include "netcdfcpp.h"
#include <map>
#include <string>

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
class CBox
{
private:
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	typedef Eigen::VectorBlock<Eigen::VectorXd,Dim> dvecBlock;

//global variables to read box configurations
	static std::map<string,CBox<Dim>*> BoxTypes;

	//enum{SQUARE, RECTANGULAR, SYMMETRIC_QUADRILATERAL};
	enum{RECTANGULAR, SYMMETRIC_QUADRILATERAL};

	int BoxSymmetry;

//Functions to check that the correct dimensions, variables, and attributes are present
//in a netCDF file.
	static bool CheckNetCDF(const NcFile &file);
	static void PopulateNetCDF(NcFile &file);
	
public:
//global functions to read box configurations
	static CBox *Read(const NcFile &file,int record);
	static void AddBoxType(string type,CBox *box);

private:
//variables specifying the transformation
	dmat Transformation; 
	dmat Inverse_Transformation;
	
public:
//constructors and copy operators
	CBox();
	CBox(const dmat Trans);
	CBox(dbl sx, dbl sy, dbl sz);
	CBox(const CBox &box);
	
	const CBox &operator=(const CBox &box);

//functions to write box configurations
	void Write(NcFile &file,int record); 
	virtual string DataToString() = 0;
	
//functions to read box configurations
 	virtual void StringToData(string Data) = 0;
 	virtual CBox *Create() = 0;
	
//functions using the transformation matrix
	void SetTransformation(dmat &Trans);
	void GetTransformation(dmat &Trans);
	void Transform(dvec &Point) const;
	void Transform(Eigen::VectorXd &Points) const;
	void InverseTransform(dvec &Point) const;
	void InverseTransform(Eigen::VectorXd &Points) const;
	void InverseTransformAndMove(Eigen::VectorXd &Points, const Eigen::VectorXd &t_Displacement);

//set and get the volume
	virtual void SetVolume(dbl V);
	virtual dbl CalculateVolume() const;

//get a list of the periodic dimensions
	virtual void GetPeriodicDimensions(std::vector<int> &) const = 0;
	
//functions involving the boundary
	virtual void MoveParticles(Eigen::VectorXd &Points,Eigen::VectorXd &Displacements)  {};
	virtual void MoveParticle(dvecBlock Point, dvec const &Displacement) {}; //dvecBlock's are themselves references, so Point should NOT be passed as a reference.
	virtual void MinimumDisplacement(const dvec &PointA, const dvec &PointB, dvec &Displacement) const {};
	virtual void MinimumDisplacement(const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointA,const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointB, dvec &Displacement) const {};

};

template <int Dim>
std::map<string,CBox<Dim>*> CBox<Dim>::BoxTypes;


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//global functions to read box configurations
template <int Dim>
CBox<Dim> *CBox<Dim>::Read(const NcFile &file,int record)
{
	if(!CheckNetCDF(file))
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a file that has no appropriate box data."));
	
	if(Dim!=file.get_dim("System_Dimension")->size())
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a record has an inconsistent nubmer of dimensions."));
	
	if(record>=file.get_dim("Records")->size())
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a record that does not exist."));

	//read the matrix
	dmat Transformation;
	NcVar *TransVar = file.get_var("Box_Transformation");
	TransVar->set_cur(record);
	TransVar->get(Transformation.data(),1,Dim,Dim);
	
	//read the box type
	char data[256];
	NcVar *DataVar = file.get_var("Box_Data");
	DataVar->set_cur(record);
	DataVar->get(data,1,256);
	string s_data = data;
	//split up the string using colons.
	vector<string> split = SplitString(s_data,":");
	//the first element of the split string should be the box name so create that kind of box
	CBox *box = BoxTypes[split[0]]->Create();
	//go set the data using the rest of the string and set the transformation
	box->StringToData(s_data);
	box->SetTransformation(Transformation);
	
	return box;
}

template <int Dim>
void CBox<Dim>::AddBoxType(string type,CBox<Dim> *box)
{
	BoxTypes[type] = box;
}

//Functions to check that the correct dimensions, variables, and attributes are present
//in a netCDF file.
template <int Dim>
bool CBox<Dim>::CheckNetCDF(const NcFile &file)
{
	return (file.get_att("Box_Populated")!=NULL);
}

//NOTE: at the moment this assumes that the box class will only ever be saved through the system class. Is this a good assumption?
template <int Dim>
void CBox<Dim>::PopulateNetCDF(NcFile &file)
{
	NcDim *recorddim = file.get_dim("Records");
	NcDim *dimdim = file.get_dim("System_Dimension");
	NcDim *datadim = file.get_dim("Data_Size");
	if(datadim==NULL)
		datadim = file.add_dim("Data_Size",256);

	file.add_att("Box_Populated",1);

	file.add_var("Box_Transformation",ncDouble,recorddim,dimdim,dimdim);
	file.add_var("Box_Data",ncChar,recorddim,datadim);
}
	
//constructors and copy operators
template <int Dim>
CBox<Dim>::CBox()
{
	Transformation = dmat::Identity();
	Inverse_Transformation = dmat::Identity();
	BoxSymmetry = SYMMETRIC_QUADRILATERAL;
}

template <int Dim>
CBox<Dim>::CBox(const dmat Trans)
{
	Transformation = Trans;
	Inverse_Transformation = Transformation.inverse();
	BoxSymmetry = SYMMETRIC_QUADRILATERAL;
}
	
template<int Dim>
CBox<Dim>::CBox(dbl sx, dbl sy, dbl sz)
{
	Transformation = dmat::Identity();
	Inverse_Transformation = Transformation.inverse();
	BoxSymmetry = SYMMETRIC_QUADRILATERAL;
}

template <int Dim>
CBox<Dim>::CBox(const CBox &box) : Transformation(box.Transformation) 
{
	Inverse_Transformation = Transformation.inverse();
	BoxSymmetry = box.BoxSymmetry;
}

template <int Dim> 
const CBox<Dim> &CBox<Dim>::operator=(const CBox<Dim> &box)
{
	Transformation = box.GetTransformation();
	Inverse_Transformation = Transformation.inverse();
	BoxSymmetry = box.BoxSymmetry;
}

//functions to write box configurations
template <int Dim>
void CBox<Dim>::Write(NcFile &file,int record)
{
	cout << "Saving box.\n";
	if(!CheckNetCDF(file))
		PopulateNetCDF(file);
		
	if(Dim!=file.get_dim("System_Dimension")->size())
		throw(CException("CBox<Dim>::WriteBox","Attempting to read a box from a record has an inconsistent nubmer of dimensions."));
	
	if(record>file.get_dim("Records")->size())
		throw(CException("CBox<Dim>::WriteBox","Attempting to read a box from a record that does not exist."));
	
	NcVar *TransVar = file.get_var("Box_Transformation");
	
	cout << "Writing transformation.\n";
	
	TransVar->set_cur(record);
	TransVar->put(Transformation.data(),1,Dim,Dim);
	
	cout << "Writing box data.\n";
	
	NcVar *DataVar = file.get_var("Box_Data");
	string str_data = DataToString();
	DataVar->set_cur(record);
	DataVar->put(str_data.c_str(),1,str_data.size());
}

	
//functions to read box configurations
	
//functions using the transformation matrix
template <int Dim>
void CBox<Dim>::SetTransformation(dmat &Trans)
{
	Transformation = Trans;
	Inverse_Transformation = Transformation.inverse();
}

template <int Dim>
void CBox<Dim>::GetTransformation(dmat &Trans)
{
	Trans = Transformation;
}

template <int Dim>
void CBox<Dim>::Transform(dvec &Point) const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			Point = Transformation.diagonal().cwiseProduct(Point); 
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			Point = Transformation * Point;	
	}
}

template <int Dim>
void CBox<Dim>::Transform(Eigen::VectorXd &Points) const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			for(int i = 0 ; i<Points.rows()/Dim ; i++)
				Points.segment<Dim>(Dim*i) = Transformation.diagonal().cwiseProduct(Points.segment<Dim>(Dim*i));
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			for(int i = 0 ; i<Points.rows()/Dim ; i++)
				Points.segment<Dim>(Dim*i) = Transformation*Points.segment<Dim>(Dim*i);
	}
}

template <int Dim>
void CBox<Dim>::InverseTransform(dvec &Point) const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			Point = Inverse_Transformation.diagonal().cwiseProduct(Point); 
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			Point = Inverse_Transformation * Point;
	}
}

template <int Dim>
void CBox<Dim>::InverseTransform(Eigen::VectorXd &Points) const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			for(int i = 0 ; i<Points.rows()/Dim ; i++)
				Points.segment<Dim>(Dim*i) = Inverse_Transformation.diagonal().cwiseProduct(Points.segment<Dim>(Dim*i));
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			for(int i = 0 ; i<Points.rows()/Dim ; i++)
				Points.segment<Dim>(Dim*i) = Inverse_Transformation*Points.segment<Dim>(Dim*i);
	}
}

template <int Dim>
void CBox<Dim>::InverseTransformAndMove(Eigen::VectorXd &Points, const Eigen::VectorXd &t_Displacement)
{
	dvec Displacement;
	int np = Points.cols()/Dim;
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			for(int i=0; i<np; ++i)
			{
				Displacement = Inverse_Transformation.diagonal().cwiseProduct(t_Displacement.segment<Dim>(Dim*i));
				MoveParticle(Points.segment<Dim>(Dim*i), Displacement);
			}
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			for(int i=0; i<np; ++i)
			{
				Displacement = Inverse_Transformation * t_Displacement.segment<Dim>(Dim*i);
				MoveParticle(Points.segment<Dim>(Dim*i), Displacement);
			}
	}
};



template <int Dim>
void CBox<Dim>::SetVolume(dbl V)
{
	dbl Vold = CalculateVolume();
	dbl Lrescale = std::pow(V/Vold,1./((dbl)Dim));
	Transformation *= Lrescale;
	Inverse_Transformation = Transformation.inverse(); //Could just divide Inverse_Transformation by Lrescale, but this is probably more stable.
}


template <int Dim>
dbl CBox<Dim>::CalculateVolume() const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			return Transformation.diagonal().prod();
		case SYMMETRIC_QUADRILATERAL:
		default:
			return fabs(Transformation.determinant());
	}
}

}

#endif
