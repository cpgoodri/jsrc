#ifndef ORDER_PARAM_H
#define ORDER_PARAM_H

#include "std_include.h"
//#include "std_fns.h"
#include "ylm.h"

class LocalOrder
{
	typedef std::complex<dbl> cdbl;
	typedef std::vector< std::vector< std::pair<int,int> > > NBR_LIST;

public:
	//Functions for 2d
	static void calc_psili_2d(std::vector< sphere_bond<2> > const &blist, NBR_LIST const &nbrs, std::vector<cdbl> &psili, int l);
	static void Calculate_Psili(std::vector< sphere_bond<2> > const &blist, int NP, std::vector<dbl> &Psili, int l=6);
	static void calc_psili_2d_normalized(std::vector< sphere_bond<2> > const &blist, NBR_LIST const &nbrs, std::vector<cdbl> &psili, int l);
	static void calc_Slij_2d(std::vector< sphere_bond<2> > const &blist, std::vector<cdbl> &psili, int l, std::vector<dbl> &Slij);

	//Functions for 3d
	static void calc_qlmi_3d(std::vector< sphere_bond<3> > const &blist, NBR_LIST const &nbrs, std::vector<cdbl> &qlmi, int l);
	static void calc_Slij_3d(std::vector< sphere_bond<3> > const &blist, std::vector<cdbl> &qlmi, int l, std::vector<dbl> &Slij);

	//Functions for both 2d and 3d
						static void Calculate_fi_from_Slij(NBR_LIST const &nbrs, std::vector<dbl> const &Slij, std::vector<dbl> &fi, dbl thresh);
	template <int Dim>	static void CalculateNbrsList(std::vector< sphere_bond<Dim> > const &blist, NBR_LIST &nbrs, int NP);
	template <int Dim>	static void CalculateSlij(std::vector< sphere_bond<Dim> > const &blist, NBR_LIST const &nbrs, int l, std::vector<dbl> &Slij);
	
	template <int Dim>	static void Calculate_fi(std::vector< sphere_bond<Dim> > const &blist, int NP, std::vector<dbl> &Slij, std::vector<dbl> &fi, dbl &fbar, dbl thresh=0.7, int l=6);
	template <int Dim>	static void Calculate_fi(std::vector< sphere_bond<Dim> > const &blist, int NP, std::vector<dbl> &fi, dbl &fbar, dbl thresh=0.7, int l=6);
	template <int Dim>	static void Calculate_fi(std::vector< sphere_bond<Dim> > const &blist, int NP, dbl &fbar, dbl thresh=0.7, int l=6);

};

////////////////////
//Functions for 2d//
////////////////////
void LocalOrder::calc_psili_2d(std::vector< sphere_bond<2> > const &blist, NBR_LIST const &nbrs, std::vector<cdbl> &psili, int l)
{
	assert(l%2==0);
	dbl theta;
	int nn, b;
	int NP = (int)nbrs.size();

	psili.assign(NP, cdbl(0.));
	for(int i=0; i<NP; ++i)
	{
		if((int)nbrs[i].size() > 0)
		{
			for(nn=0; nn<(int)nbrs[i].size(); ++nn)
			{
				b = nbrs[i][nn].second;
				theta = atan2(blist[b].r[1], blist[b].r[0]);
				printf("theta = % f\n", theta/M_PI);
				psili[i] += cdbl(std::cos(l*theta),std::sin(l*theta));
			}
			psili[i] /= (dbl)nbrs[i].size();
		}
	}
}

void LocalOrder::Calculate_Psili(std::vector< sphere_bond<2> > const &blist, int NP, std::vector<dbl> &Psili, int l)
{
	//Set nbrs list
	NBR_LIST nbrs;
	CalculateNbrsList(blist, nbrs, NP);

	std::vector<cdbl> psili;
	calc_psili_2d(blist, nbrs, psili, l);

	Psili.assign(NP, 0.);
	assert(Psili.size() == psili.size());
	for(int i=0; i<NP; ++i)
		Psili[i] = std::norm(psili[i]);
}

void LocalOrder::calc_psili_2d_normalized(std::vector< sphere_bond<2> > const &blist, NBR_LIST const &nbrs, std::vector<cdbl> &psili, int l)
{
	assert(l%2==0);
	dbl theta;
	int nn, b;
	int NP = (int)nbrs.size();
	dbl mag,mag2;

	psili.assign(NP, cdbl(0.));
	for(int i=0; i<NP; ++i)
	{
		if((int)nbrs[i].size() > 0)
		{
			for(nn=0; nn<(int)nbrs[i].size(); ++nn)
			{
				b = nbrs[i][nn].second;
				theta = atan2(blist[b].r[1], blist[b].r[0]);
				psili[i] += cdbl(std::cos(l*theta),std::sin(l*theta));
			}
			//normalize
			mag = std::sqrt( std::norm(psili[i]) );
			psili[i] /= mag;
		}
	}
}

void LocalOrder::calc_Slij_2d(std::vector< sphere_bond<2> > const &blist, std::vector<cdbl> &psili, int l, std::vector<dbl> &Slij)
{
	Slij.assign((int)blist.size(), 0.);

	int i, j;
	for(int b=0; b<(int)blist.size(); ++b)
	{
		i = blist[b].i;
		j = blist[b].j;
		Slij[b] += std::real(psili[i] * std::conj(psili[j]) );
	}
}


template <>
void LocalOrder::CalculateSlij<2>(std::vector< sphere_bond<2> > const &blist, NBR_LIST const &nbrs, int l, std::vector<dbl> &Slij)
{
	int NP = (int)nbrs.size();

	//Calculate the psi_l's
	std::vector<cdbl> psili;
	calc_psili_2d_normalized(blist, nbrs, psili, l);
	assert((int)psili.size() == NP);

	//Calculate the Slij's
	calc_Slij_2d(blist, psili, l, Slij);
}


////////////////////
//Functions for 3d//
////////////////////
void LocalOrder::calc_qlmi_3d(std::vector< sphere_bond<3> > const &blist, NBR_LIST const &nbrs, std::vector<cdbl> &qlmi, int l)
{
	const int c2lp1 = 2*l+1;
	assert(l%2==0);

	dbl theta, phi;
	int nn,m,b;
	int NP = (int)nbrs.size();
	dbl mag,mag2;

	qlmi.assign(NP*c2lp1, cdbl(0.));
	for(int i=0; i<NP; ++i)
	{
		if((int)nbrs[i].size() > 0)
		{
			for(nn=0; nn<(int)nbrs[i].size(); ++nn)
			{
				b = nbrs[i][nn].second;
				theta = acos(blist[b].r[2]/blist[b].drlen);
				phi = atan2(blist[b].r[1], blist[b].r[0]);
				for(m=-l; m<=l; ++m)
					qlmi[i*c2lp1 + m+l] += Ylm::calc_Ylm(l,m, theta, phi);
			}
			//normalize
			mag2 = 0.;
			for(m=-l; m<=l; ++m)
				mag2 += std::norm(qlmi[i*c2lp1 + m+l]);
			mag = std::sqrt(mag2);
			for(m=-l; m<=l; ++m)
				qlmi[i*c2lp1 + m+l] /= mag;
		}
	}
}

void LocalOrder::calc_Slij_3d(std::vector< sphere_bond<3> > const &blist, std::vector<cdbl> &qlmi, int l, std::vector<dbl> &Slij)
{
	const int c2lp1 = 2*l+1;
	Slij.assign((int)blist.size(), 0.);

	int i, j, m;
	for(int b=0; b<(int)blist.size(); ++b)
	{
		i = blist[b].i;
		j = blist[b].j;
		for(m=-l; m<=l; ++m)
			Slij[b] += std::real(qlmi[i*c2lp1 + m+l] * std::conj(qlmi[j*c2lp1 + m+l]));
		//Slij[b] /= calc_qlmi_mag(qlmi, l, i) * calc_qlmi_mag(qlmi, l, j); //Dont need this because qlmi is already normalized.
	}
}

template <>
void LocalOrder::CalculateSlij<3>(std::vector< sphere_bond<3> > const &blist, NBR_LIST const &nbrs, int l, std::vector<dbl> &Slij)
{
	int NP = (int)nbrs.size();

	//Calculate the q_lm's
	std::vector<cdbl> qlmi;
	calc_qlmi_3d(blist, nbrs, qlmi, l);
	assert((int)qlmi.size() == (2*l+1)*NP);

	//Calculate the Slij's
	calc_Slij_3d(blist, qlmi, l, Slij);
}





////////////////////////////////
//Functions for both 2d and 3d//
////////////////////////////////

template <int Dim>
void LocalOrder::CalculateNbrsList(std::vector< sphere_bond<Dim> > const &blist, NBR_LIST &nbrs, int NP)
{
	nbrs.assign(NP, std::vector< std::pair<int,int> >());
	int i,j, Nc = (int)blist.size();
	for(int b=0; b<(int)Nc; ++b)
	{
		i = (int)blist[b].i;
		j = (int)blist[b].j;
		nbrs[i].push_back(std::make_pair(j,b));
		nbrs[j].push_back(std::make_pair(i,b));
	}
	assert((int)nbrs.size()==NP);
}

void LocalOrder::Calculate_fi_from_Slij(NBR_LIST const &nbrs, std::vector<dbl> const &Slij, std::vector<dbl> &fi, dbl thresh)
{
	int num_connected;
	int NP = (int)nbrs.size();
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

template <int Dim>
void LocalOrder::Calculate_fi(std::vector< sphere_bond<Dim> > const &blist, int NP, std::vector<dbl> &Slij, std::vector<dbl> &fi, dbl &fbar, dbl thresh, int l)
{
	int Nc = (int)blist.size();

	//Set nbrs list
	NBR_LIST nbrs;
	CalculateNbrsList(blist, nbrs, NP);

	//Calculate the Slij's
	CalculateSlij(blist, nbrs, l, Slij);
	assert((int)Slij.size() == Nc);

	Calculate_fi_from_Slij(nbrs, Slij, fi, thresh);
	
	calc_avg(fi, fbar);
}

template <int Dim>
void LocalOrder::Calculate_fi(std::vector< sphere_bond<Dim> > const &blist, int NP, std::vector<dbl> &fi, dbl &fbar, dbl thresh, int l)
{
	std::vector<dbl> Slij;
	Calculate_fi(blist, NP, Slij, fi, fbar, thresh, l);
}

template <int Dim>
void LocalOrder::Calculate_fi(std::vector< sphere_bond<Dim> > const &blist, int NP, dbl &fbar, dbl thresh, int l)
{
	std::vector<dbl> Slij;
	std::vector<dbl> fi;
	Calculate_fi(blist, NP, Slij, fi, fbar, thresh, l);
}








#endif //ORDER_PARAM_H
