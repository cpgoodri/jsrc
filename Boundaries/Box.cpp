//global functions to read box configurations
template <int Dim>
static CBox<Dim> *CBox<Dim>::ReadBox(const NcFile &file,int record)
{
	if(!CheckNetCDF(file))
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a file that has no appropriate box data."));
	
	if(Dim!=file.get_dim("Dimension")->size())
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a record has an inconsistent nubmer of dimensions."));
	
	if(record>=file.get_dim("Records")->size())
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a record that does not exist."));

	//read the matrix
	Eigen::Matrix<double,Dim,Dim> Transformation;
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
static void CBox<Dim>::AddBoxType(string type,CBox<Dim> *box)
{
	BoxTypes.insert(type,box);
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

	file.add_var("Box_Transformation",recorddim,dimdim,dimdim);
	file.add_var("Box_Data",recorddim,datadim);
}
	
//constructors and copy operators
template <int Dim>
CBox<Dim>::CBox()
{
	Transformation = Matrix<double,Dim,Dim>::Identity();
	Inverse_Transformation = Matrix<double,Dim,Dim>::Identity();
}

template <int Dim>
CBox<Dim>::CBox(const Eigen::Matrix<double,Dim,Dim> Trans)
{
	Transformation = Trans;
	Inverse_Transformation = Transformation.inverse();
}
	
template<int Dim>
CBox<Dim>::CBox(double sx, double sy, double sz)
{
	Transformation = Matrix<double,Dim,Dim>::Identity();
	Inverse_Transformation = Transformation.inverse();
}

template <int Dim>
CBox<Dim>::CBox(const CBox &box) : Transformation(box.Transformation) 
{
	Inverse_Transformation = Transformation.inverse();
}

template <int Dim> 
const CBox &CBox<Dim>::operator=(const CBox &box)
{
	Transformation = box.GetTransformation();
	Inverse_Transformation = Transformation.inverse();
}

//functions to write box configurations
template <int Dim>
void CBox<Dim>::WriteBox(NcFile &file,int record)
{
	if(!CheckNetCDF(file))
		PopulateNetCDF(file);
		
	if(Dim!=file.get_dim("Dimension")->size())
		throw(CException("CBox<Dim>::WriteBox","Attempting to read a box from a record has an inconsistent nubmer of dimensions."));
	
	if(record>file.get_dim("Records")->size())
		throw(CException("CBox<Dim>::WriteBox","Attempting to read a box from a record that does not exist."));
	
	NcVar *TransVar = file.get_var("Box_Transformation");
	
	TransVar->set_cur(record);
	TransVar->put(Transformation.data(),1,Dim,Dim);
	
	NcVar *DataVar = file.get_var("Box_Data");
	string str_data = DataToString();
	DataVar->set_cur(record);
	DataVar->put(str_data.c_str(),1,str_data.size());
}

	
//functions to read box configurations
	
//functions using the transformation matrix
template <int Dim>
void CBox<Dim>::SetTransformation(Eigen::Matrix<double,Dim,Dim> &Trans)
{
	Transformation = Trans;
	Inverse_Transformation = Transformation.inverse();
}

template <int Dim>
void CBox<Dim>::GetTransformation(Eigen::Matrix<double,Dim,Dim> &Trans)
{
	return Transformation;
}

template <int Dim>
void CBox<Dim>::Transform(Eigen::Matrix<double,Dim,1> &Point)
{
	Point = Transformation * Point;	
}

template <int Dim>
void CBox<Dim>::Transform(Eigen::VectorXd &Points)
{
	for(int i = 0 ; i<Points.cols()/Dim ; i++)
		Points.segment<Dim>(Dim*i) = Transformation*Points.segment<Dim>(Dim*i);
}

template <int Dim>
void CBox<Dim>::InverseTransform(Eigen::Matrix<double,Dim,1> &Point)
{
	Point = Inverse_Transformation * Point;
}

template <int Dim>
void CBox<Dim>::InverseTransform(Eigen::VectorXd &Points)
{
	for(int i = 0 ; i<Points.cols()/Dim ; i++)
		Points.segment<Dim>(Dim*i) = Inverse_Transformation*Points.segment<Dim>(Dim*i);
}