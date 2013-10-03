#ifndef UMFPACK_INTERFACE_H
#define UMFPACK_INTERFACE_H

#ifndef DONT_USE_SUITESPARSE

#include "../Resources/std_include.h"
#include "umfpack.h"
#include <Eigen/Sparse>
#include "CCSmatrix.h"


#define UMF_USE_LONG

#ifdef UMF_USE_LONG
	typedef long UMFint;
#else
	typedef int UMFint;
#endif


class UmfpackInterface
{
public:
	//Symbolic and Numeric are void* objects used by UMFPACK
	void *Symbolic;
	void *Numeric;

	//Eigen::SparseMatrix<dbl> const *A;
	CCSmatrix<UMFint> A;

	//Allocation flags.
	bool symbolic_allocated;
	bool numeric_allocated;
	bool symbolic_found;
	bool numeric_found;

	//options
	int  verbose_invert;		//default 0

private:
	//other internal data
	
public:
	//Constructors/destructors
	UmfpackInterface();
	~UmfpackInterface();

private:
	void clear_symbolic();
	void clear_numeric();
	void init();

public:
	void set_default_options();
	void set_verbose_invert     (int yesno ){	verbose_invert      = yesno;	};

	int info() const;
	int compute(Eigen::SparseMatrix<dbl> const &_A);
	Eigen::VectorXd solve(Eigen::VectorXd const &B);

	UMFint UMF_symbolic(UMFint n_row, UMFint n_col, const UMFint Ap[], const UMFint Ai[], const double Ax[], void **Symbolic, const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO]);
	UMFint UMF_numeric (const UMFint Ap [ ], const UMFint Ai [ ], const double Ax [ ], void *Symbolic, void **Numeric, const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]);
	UMFint UMF_solve (UMFint sys, const UMFint Ap [ ], const UMFint Ai [ ], const double Ax [ ], double X [ ], const double B [ ], void *Numeric, const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]);
	void UMF_free_symbolic(void **Symbolic);
	void UMF_free_numeric (void **Numeric );
	void UMF_defaults (double Control [UMFPACK_CONTROL]);
};

#endif //DONT_USE_SUITESPARSE


#endif //UMFPACK_INTERFACE_H

