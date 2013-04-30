#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "../Resources/std_include.h"
#include "../Computers/BaseComputer.h"
#include "minimizerFIRE.h"

namespace LiuJamming
{

template<int Dim>
class CSimpleMinimizer
{
private:
	enum{NONE=0,FIRE, LBFGS, CG_PR, CG_FR};
	CBaseComputer<Dim> *pComputer;
	int N_dof;

public:
	CSimpleMinimizer(CBaseComputer<Dim> &TargetComputer, int minimize)
		: pComputer(&TargetComputer)
	{
		N_dof = pComputer->GetNdof();
		switch(minize)
			minimizeFIRE();
	};

	CSimpleMinimizer(CBaseComputer<Dim> &TargetComputer, dbl tol=1e-12, dbl delta_t_start=0.1, int max_iterations=-1)
		: pComputer(&TargetComputer)
	{
		N_dof = pComputer->GetNdof();
		minimizeFIRE(
	};

	void minimizeFIRE(dbl tol, dbl delta_t_start, int max_iterations = -1);
};

template<int Dim>
void CSimpleMinimizer<Dim>::minimizeFIRE(dbl tol, dbl delta_t_start, int max_iterations)
{
	typedef CMinimizerFIRE< CBaseComputer<Dim> > MINIMIZER;
	MINIMIZER mm(N_dof);

	mm.set_functions(pComputer, &CBaseComputer<Dim>::Evaluate, &CBaseComputer<Dim>::Progress, &CBaseComputer<Dim>::Move, &CBaseComputer<Dim>::ReportHeader);

	if(max_iterations > 0)
		mm.set_max_iterations(max_iterations);

//	dbl delta_t_start = 0.04*plist.get_length_unit();
	mm.set_FIRE_delta_t_start(delta_t_start);
	mm.set_FIRE_delta_t_max(10.*delta_t_start);

	mm.minimize(tol);
};


}

#endif //MINIMIZER
