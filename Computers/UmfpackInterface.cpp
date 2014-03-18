
#ifndef DONT_USE_SUITESPARSE

#include "UmfpackInterface.h"

UmfpackInterface::UmfpackInterface()
{
	init();
}

UmfpackInterface::~UmfpackInterface()
{
	clear_symbolic();
	clear_numeric();
}

void UmfpackInterface::clear_symbolic()
{
	if(!symbolic_allocated) return;
	UMF_free_symbolic(&Symbolic); //This might change for different matrix types.
	symbolic_allocated = false;
}

void UmfpackInterface::clear_numeric()
{
	if(!numeric_allocated) return;
	UMF_free_numeric(&Numeric); //This might change for different matrix types.
	numeric_allocated = false;
}

void UmfpackInterface::init()
{
	symbolic_allocated = false;
	numeric_allocated = false;
	symbolic_found = false;
	numeric_found = false;

	set_default_options();
}

void UmfpackInterface::set_default_options()
{
	verbose_invert = 0;
}


int UmfpackInterface::info() const
{
	if(symbolic_found && numeric_found)
		return Eigen::Success;
	return Eigen::NumericalIssue;
}

int UmfpackInterface::compute(Eigen::SparseMatrix<dbl> const &_A)
{
	A.SetMatrix(_A);
	UMFint const *Ap = A.Ap.data();
	UMFint const *Ai = A.Ai.data();
	dbl    const *Ax = A.Ax.data();

	dbl Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
	UMF_defaults(Control);
	Control[UMFPACK_PIVOT_TOLERANCE] = 1.0; //Use true partial pivoting.
	//Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;

    //Calcuate the symbolic decomposition.
	if(UMF_symbolic(A.n_rows, A.n_cols, Ap, Ai, Ax, &Symbolic, Control, Info) != UMFPACK_OK)	printf("Error, symbolic failed\n");


	switch ((UMFint)Info[UMFPACK_STATUS])
	{
		case UMFPACK_OK:if(verbose_invert)			printf("UMFPACK symbolic factorization completed successfully in %e seconds", Info[UMFPACK_SYMBOLIC_TIME]); break;
		case UMFPACK_ERROR_invalid_matrix:			printf("Error: invalid matrix\n"); break;
		case UMFPACK_ERROR_out_of_memory:			printf("Error: UMFPACK out of memory\n"); break;
		case UMFPACK_ERROR_invalid_Symbolic_object:	printf("Error: invalid Symbolic object\n"); break;
		case UMFPACK_ERROR_internal_error:			printf("Serious Error: Please contact UMFPACK author.\n"); break;
		default:									printf("Unknown Error.\n"); break;
	}
	if(verbose_invert)switch ((UMFint)Info[UMFPACK_STRATEGY_USED])
	{
		case UMFPACK_STRATEGY_SYMMETRIC:			printf("\tStrategy: symmetric\n"); break;
		case UMFPACK_STRATEGY_UNSYMMETRIC:			printf("\tStryatey: unsymmetric\n"); break;
		default:									printf("\tError: unknown symbolic strategy...\n"); break;
	}
	symbolic_allocated = 1;
	if(Info[UMFPACK_STATUS]==UMFPACK_OK)
		symbolic_found = 1;
	else
		symbolic_found = 0;

	//Calculate the numerical decomposition.
	if(UMF_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info) != UMFPACK_OK)		printf("Error, numeric failed\n");
	switch ((UMFint)Info[UMFPACK_STATUS])
	{
		case UMFPACK_OK: if(verbose_invert) 		printf("UMFPACK numeric inversion completed successfully in %e seconds", Info[UMFPACK_NUMERIC_TIME]); break;
		case UMFPACK_WARNING_singular_matrix:		printf("Warning: singular matrix\n"); break;
		case UMFPACK_ERROR_out_of_memory:			printf("Error: UMFPACK out of memory\n"); break;
		case UMFPACK_ERROR_invalid_Symbolic_object:	printf("Error: invalid Symbolic object\n"); break;
		default:									printf("Unknown Error.\n"); break;
	}
	if(verbose_invert)								printf("\tComputation used %e FLOPS and had a peak memory usage of %e\n", Info[UMFPACK_FLOPS], Info[UMFPACK_PEAK_MEMORY]);
	numeric_allocated = 1;
	if(Info[UMFPACK_STATUS]==UMFPACK_OK)
		numeric_found = 1;
	else
		numeric_found = 0;

	return info();
}




Eigen::VectorXd UmfpackInterface::solve(Eigen::VectorXd const &B)
{
	if(!numeric_found)
	{
		printf("Warning, trying to solve Mx=b without a valid numeric object.\n"); fflush(stdout);
		exit(EXIT_FAILURE);
	}

	UMFint const *Ap = A.Ap.data();
	UMFint const *Ai = A.Ai.data();
	dbl    const *Ax = A.Ax.data();
	
	Eigen::VectorXd X = Eigen::VectorXd::Zero(B.size());
	dbl *x = X.data();
	dbl const *b = B.data();

	double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
	UMF_defaults(Control);
	Control[UMFPACK_IRSTEP] = 1000;

	if(UMF_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, Control, Info) != UMFPACK_OK)		printf("Error, solve failed\n");

	switch ((UMFint)Info[UMFPACK_STATUS])
	{
		case UMFPACK_OK: if(verbose_invert)			printf("UMFPACK system solved successfully in %e seconds\t\t", Info[UMFPACK_SOLVE_TIME]); break;
		case UMFPACK_WARNING_singular_matrix:		printf("Warning: singular matrix\n"); break;
		case UMFPACK_ERROR_out_of_memory:			printf("Error: UMFPACK out of memory\n"); break;
		case UMFPACK_ERROR_invalid_Numeric_object:	printf("Error: invalid Numeric object\n"); break;
		default:									printf("Unknown Error.\n"); break;
	}
	if(verbose_invert)								printf("\tComputation used %e FLOPS\n", Info[UMFPACK_SOLVE_FLOPS]);


	return X;
//	return (Info[UMFPACK_STATUS]==UMFPACK_OK)?1:0;
}

#ifdef UMF_USE_LONG

	UMFint UmfpackInterface::UMF_symbolic(UMFint n_row, UMFint n_col, const UMFint Ap[], const UMFint Ai[], const double Ax[], void **Symbolic, const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO])
	{
		return umfpack_dl_symbolic(n_row, n_col, Ap, Ai, Ax, Symbolic, Control, Info);
	}
		
	UMFint UmfpackInterface::UMF_numeric (const UMFint Ap [ ], const UMFint Ai [ ], const double Ax [ ], void *Symbolic, void **Numeric, const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
	{
		return umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, Numeric, Control, Info);
	}

	UMFint UmfpackInterface::UMF_solve (UMFint sys, const UMFint Ap [ ], const UMFint Ai [ ], const double Ax [ ], double X [ ], const double B [ ], void *Numeric, const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
	{
		return umfpack_dl_solve(sys, Ap, Ai, Ax, X, B, Numeric, Control, Info);
	}

	void UmfpackInterface::UMF_free_symbolic(void **Symbolic)
	{
		umfpack_dl_free_symbolic(Symbolic);
	}

	void UmfpackInterface::UMF_free_numeric (void **Numeric )
	{
		umfpack_dl_free_numeric(Numeric);
	}

	void UmfpackInterface::UMF_defaults (double Control [UMFPACK_CONTROL])
	{
		umfpack_dl_defaults(Control);
	}

#else

	UMFint UmfpackInterface::UMF_symbolic(UMFint n_row, UMFint n_col, const UMFint Ap[], const UMFint Ai[], const double Ax[], void **Symbolic, const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO])
	{
		return umfpack_di_symbolic(n_row, n_col, Ap, Ai, Ax, Symbolic, Control, Info);
	}
		
	UMFint UmfpackInterface::UMF_numeric (const UMFint Ap [ ], const UMFint Ai [ ], const double Ax [ ], void *Symbolic, void **Numeric, const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
	{
		return umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, Control, Info);
	}

	UMFint UmfpackInterface::UMF_solve (UMFint sys, const UMFint Ap [ ], const UMFint Ai [ ], const double Ax [ ], double X [ ], const double B [ ], void *Numeric, const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
	{
		return umfpack_di_solve(sys, Ap, Ai, Ax, X, B, Numeric, Control, Info);
	}

	void UmfpackInterface::UMF_free_symbolic(void **Symbolic)
	{
		umfpack_di_free_symbolic(Symbolic);
	}

	void UmfpackInterface::UMF_free_numeric (void **Numeric )
	{
		umfpack_di_free_numeric(Numeric);
	}

	void UmfpackInterface::UMF_defaults (double Control [UMFPACK_CONTROL])
	{
		umfpack_di_defaults(Control);
	}

#endif 



#endif //DONT_USE_SUITESPARSE


