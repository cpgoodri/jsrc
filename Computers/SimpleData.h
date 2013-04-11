#ifndef SIMPLE_DATA_H
#define SIMPLE_DATA_H

template<int Dim>
class CSimpleData
{
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
public:
	int NPp;
	int Nc;
	int NcmNciso;
	dbl Volume;
	dbl energy;
	dbl pressure;
	dmat stress;
	dbl max_grad;

	inline dbl Z()      const { return 2.*((dbl)Nc)/((dbl)NPp); };
	inline dbl deltaZ() const { return 2.*((dbl)NcmNciso)/((dbl)NPp); };
};

#endif //SIMPLE_DATA_H

