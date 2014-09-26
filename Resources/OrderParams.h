#ifndef ORDER_PARAM_H
#define ORDER_PARAM_H

#include "std_include.h"
#include "Exception.h"
//#include "std_fns.h"
#include "Ylm.h"
//#include "../Computers/BondList.h"


namespace LiuJamming
{

template <int Dim> class CBondList;

template <int Dim>
class CLocalOrder
{
	typedef vector< vector< std::pair<int,int> > > NBR_LIST;

public:
	//Input info
	const CBondList<Dim> *pbonds;
	const int NP;
	const int Nbonds;
	const int l;
	const int num_m;

	//Calculated info
	NBR_LIST nbrs;
	vector<cdbl> qlmi;
	vector< dbl> Slij;
	vector< dbl> fi;


	CLocalOrder(CBondList<Dim> const &_bonds, int _l=6);
private:
	void init();

public:
	void Calculate_qlmi();
	void Calculate_Slij();
	void Calculate_fi(dbl thresh=0.7);
	void Calculate(dbl thresh=0.7);
};

static int set_num_m(int Dim, int _l)
{
	if(Dim==2) return 1;
	if(Dim==3) return 2*_l+1;
	throw( CException("set_num_m","Cannot construct a CLocalOrder object for the requested dimension.") );
	return 0;
}

template <int Dim> CLocalOrder<Dim>::CLocalOrder(CBondList<Dim> const &_bonds, int _l)
	: pbonds(&_bonds), NP(pbonds->GetN()), l(_l), Nbonds(pbonds->GetNBonds()), num_m(set_num_m(Dim,_l))
{
	init();
}

template <int Dim>
void CLocalOrder<Dim>::init()
{
	assert(l%2==0);

	//Set nbrs list
	nbrs.assign(NP, vector< std::pair<int,int> >());
	int i,j;
	for(int b=0; b<Nbonds; ++b)
	{
		i = (*pbonds)[b].i;
		j = (*pbonds)[b].j;
		assert(i>=0);
		assert(j>=0);
		assert(i<NP);
		assert(j<NP);

		nbrs[i].push_back(std::make_pair(j,b));
		nbrs[j].push_back(std::make_pair(i,b));
	}
	assert((int)nbrs.size()==NP);
	
	//initialize qlmi, Slij, and fi
//	qlmi.assign(NP*num_m, cdbl(0.,0.));
//	Slij.assign(Nbonds, 0.);
//	fi  .assign(NP, 0.);
}


template <> void CLocalOrder<2>::Calculate_qlmi()
{
	dbl theta;
	int nn;
	CBond<2> bnd;

	qlmi.assign(NP, cdbl(0.));
	for(int i=0; i<NP; ++i)
	{
		if((int)nbrs[i].size() > 0)
		{
			for(nn=0; nn<(int)nbrs[i].size(); ++nn)
			{
				pbonds->GetBond(nbrs[i][nn].second, bnd);
				theta = atan2(bnd.r[1], bnd.r[0]);
				qlmi[i] += cdbl(std::cos(l*theta),std::sin(l*theta));
			}
			qlmi[i] /= (dbl)nbrs[i].size();
		}
	}
}

template <> void CLocalOrder<3>::Calculate_qlmi()
{
	dbl theta, phi;
	int nn,m,b;
	dbl mag,mag2;

	CBond<3> bnd;

	qlmi.assign(NP*num_m, cdbl(0.));
	for(int i=0; i<NP; ++i)
	{
		if((int)nbrs[i].size() > 0)
		{
			for(nn=0; nn<(int)nbrs[i].size(); ++nn)
			{
				pbonds->GetBond(nbrs[i][nn].second, bnd);
				theta = acos(bnd.r[2]/bnd.rlen);
				phi = atan2(bnd.r[1], bnd.r[0]);
				for(m=-l; m<=l; ++m)
					qlmi[i*num_m + m+l] += Ylm::calc_Ylm(l,m, theta, phi);
			}
			for(m=-l; m<=l; ++m)
				qlmi[i*num_m + m+l] /= (dbl)nbrs[i].size();
		}
	}
}

static dbl normalized_product(vector<cdbl> &vec, int i, int j, int lng)
{
	Eigen::VectorXcd vec1 = Eigen::VectorXcd::Zero(lng);
	Eigen::VectorXcd vec2 = Eigen::VectorXcd::Zero(lng);

	for(int ii=0; ii<lng; ++ii)
	{
		vec1[ii] = vec[i*lng+ii];
		vec2[ii] = vec[j*lng+ii];
	}
	return std::real(vec1.dot(vec2))/(vec1.norm()*vec2.norm());
}

template<int Dim> void CLocalOrder<Dim>::Calculate_Slij()
{
	Slij.assign(Nbonds, 0.);
	int i, j, m;
	for(int b=0; b<Nbonds; ++b)
	{
		i = (*pbonds)[b].i;
		j = (*pbonds)[b].j;
		Slij[b] = normalized_product(qlmi, i, j, num_m);
	}
}

template<int Dim> void CLocalOrder<Dim>::Calculate_fi(dbl thresh)
{
	int num_connected;
	fi.assign(NP, 0.);

	for(int i=0; i<NP; ++i)
	{
		if(nbrs[i].size() > 0)
		{
			num_connected = 0;
			for(int nn=0; nn<(int)nbrs[i].size(); ++nn)
				if(Slij[nbrs[i][nn].second] > thresh)
					++num_connected;
			fi[i] = (((dbl)num_connected)/((dbl)nbrs[i].size()));
		}
	}
}


template<int Dim> void CLocalOrder<Dim>::Calculate(dbl thresh)
{
	Calculate_qlmi();
	Calculate_Slij();
	Calculate_fi(thresh);
}









}



#endif //ORDER_PARAM_H
