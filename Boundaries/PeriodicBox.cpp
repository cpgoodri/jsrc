//constructors and copy operators
template <int Dim>
CPeriodicBox<Dim>::CPeriodicBox() : CBox()
{

}

template <int Dim>
CPeriodicBox<Dim>::CPeriodicBox(const Eigen::Matrix<double,Dim,Dim> Trans) : CBox(Trans)
{


}

template <int Dim>
CPeriodicBox<Dim>::CPeriodicBox(double sx, double sy, double sz) : CBox(sx,sy,sz)
{

}

template <int Dim>
CPeriodicBox<Dim>::CPeriodicBox(const CPeriodicBox &box) : CBox(box)
{

}

template <int Dim>	
const CPeriodicBox<Dim>::CPeriodicBox &operator=(const CPeriodicBox &box)
{
	CBox::operator=(box);
	return *this;
}

//functions to write box configurations
template <int Dim>
string CPeriodicBox<Dim>::DataToString()
{
	return "PeriodicBox";
}
	
//functions to read box configurations
template <int Dim>
void CPeriodicBox<Dim>::StringToData(string Data)
{

}

template <int Dim>
CBox *CPeriodicBox<Dim>::Create()
{
	return new PeriodicBox();
}
 	
//functions involving the boundary
template <int Dim>
void CPeriodicBox<Dim>::MoveParticles(Eigen::VectorXd &Points,const Eigen::VectorXd &Displacements)
{
	Points += Displacements;
	for(int i = 0; i< Points.cols() ; i++)
		Points(i) -= floor(Points(i));
}

template <int Dim>
void CPeriodicBox<Dim>::MinimumDisplacement(const Eigen::Matrix<double,Dim,1> &PointA,const Eigen::Matrix<double,Dim,1> &PointB, Eigen::Matrix<double,Dim,1> &Displacement)
{
	Displacement = PointA-PointB;
	for(int i = 0 ; i < Dim ; i++)
		if(abs(Displacement(i))>0.5)
			Displacement(i)-=sgn(Displacement(i));
}