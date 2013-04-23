#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "../Resources/std_include.h"
#include "../Computers/BaseComputer.h"
#include "minimizerFIRE.h"

template<int Dim>
class CSimpleMinimizer
{
private:
	CBaseComputer<Dim> *pComputer;
	int N_dof;

public:
	CSimpleMinimizer(CBaseComputer<Dim> &TargetComputer, _N_dof)
		: pComputer(&TargetComputer), N_dof(_N_dof)
	{};

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

#endif //MINIMIZER
