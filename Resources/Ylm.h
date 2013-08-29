#ifndef YLM_H
#define YLM_H

#include "std_include.h"
//#include "std_fns.h"


//Calculate the associated Legendre polynomial Plm(x)
class Plm
{
public:
	static dbl calc_Plm(int l, int m, dbl x)
	{
		assert(0 <= m);
		assert(m <= l);
	
		dbl somx2, fact;
		dbl Pmm = 1.;

		//Compute Pmm
		if(m>0)
		{
			somx2 = sqrt((1.-x)*(1.+x));
			fact = 1.;
			for(int i=1; i<=m; ++i)
			{
				Pmm *= -fact*somx2;
				fact += 2.;
			}
		}
		if(l == m) return Pmm;
		
		//Compute Pm,m+1
		dbl Pmmp1 = x*((dbl)(2*m+1))*Pmm;
		if(l == m+1) return Pmmp1;

		//Compute Pm,l
		dbl Pll;
		for(int ll = m+2; ll<=l; ++ll)
		{
			Pll = (x*((dbl)(2*ll-1))*Pmmp1 - ((dbl)(ll+m-1))*Pmm)/((dbl)(ll-m));
			Pmm = Pmmp1;
			Pmmp1 = Pll;
		}
		return Pll;
	};
};

//calculate the spherical harmonic Ylm(theta, phi)
class Ylm
{
public:
	static cdbl calc_Ylm(int l, int m, dbl theta, dbl phi)
	{
		assert(l>=0);
		assert(m <= l);
		assert(-l <= m);

		if(m<0)
		{
			return (((-m)%2==1)? -1.:1.) * std::conj(calc_Ylm(l,-m,theta,phi));
		}
		assert(m>=0);

		dbl l_minus_m_factorial = factorial(l-m);
		dbl l_plus_m_factorial  = factorial(l+m);
		dbl sqrt_part = sqrt( (((dbl)(2*l+1))/(4.*M_PI)) * l_minus_m_factorial / l_plus_m_factorial );
		dbl plm = Plm::calc_Plm(l, fabs(m), cos(theta));

		return sqrt_part * plm * exp( (cdbl(0.,1.)) * ((dbl)m) * phi );
	};

	static dbl calc_Ylm_re(int l, int m, dbl theta, dbl phi)
	{
		return (calc_Ylm(l,m,theta,phi)).real();
	};

	static dbl calc_Ylm_im(int l, int m, dbl theta, dbl phi)
	{
		return (calc_Ylm(l,m,theta,phi)).imag();
	};
};







































#endif //YLM_H
