#ifndef MATRIX_INTERFACE_H
#define MATRIX_INTERFACE_H

#include "../Resources/std_include.h"

#include "arch.h"
#include "arscomp.h"
#include "arssym.h"

#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>


namespace LiuJamming
{



template <typename T> class MatrixInterface;

//Functions to interface with ARPACK++
template <typename T> class ARWrapper {}; 

template <>
class ARWrapper<dbl> 
{
public: 
	ARSymStdEig<dbl,MatrixInterface<dbl> > Solver;
	static const std::string which;
	static const std::string which2;
};
const std::string ARWrapper<dbl>::which("SA");
const std::string ARWrapper<dbl>::which2("LA");

template <>
class ARWrapper<cdbl> 
{
public: 
	ARCompStdEig<cdbl,MatrixInterface<cdbl> > Solver;
	static const std::string which;
	static const std::string which2;
};
const std::string ARWrapper<cdbl>::which("SM");
const std::string ARWrapper<cdbl>::which2("LM");



template <typename T>
static void set_identity(Eigen::SparseMatrix<T> &Identity, int nrows)
{
	assert(Identity.rows() == nrows);
	assert(Identity.cols() == nrows);

	std::vector< Eigen::Triplet<T> > triplets;
	triplets.reserve(nrows);
	for(int i=0; i<nrows; ++i)
		triplets.push_back(Eigen::Triplet<T>(i,i,((T)1.)));
	Identity.setFromTriplets(triplets.begin(), triplets.end());
}







template<typename T>
class MatrixInterface 
{
private:
	typedef Eigen::SparseMatrix<T> EMatrix;

	//internal data
	double diagonalization_time;
	long int Mv_counter;
	long int num_Mv_calls;

	long int OPv_counter;
	long int num_OPv_calls;

public:
	EMatrix A; //The matrix

//	//Used for diagonalizing with shift and invert mode
//	EMatrix OPinv; //OP = (A-sigma.I)^{-1}, where sigma is a shift and I is the identity. OPinv = A-sigma.I will be allocated only when needed.
//	Eigen::UmfPackLU<EMatrix> OP_solver;

	T *Eigenvalues;
	T *Eigenvectors;
	int num_request;
	int num_converged;
	bool compute_vecs;
	bool verbose_diagonalize;

	//Solver for UMFPACK
	Eigen::UmfPackLU<EMatrix> UMFsolver;

	MatrixInterface()
		: num_request(0),num_converged(0),compute_vecs(1),verbose_diagonalize(1), Eigenvalues(NULL), Eigenvectors(NULL)
	{
	};

	MatrixInterface(EMatrix const &mat)
		: num_request(0),num_converged(0),compute_vecs(1),verbose_diagonalize(1), Eigenvalues(NULL), Eigenvectors(NULL)
		  A(mat)
	{
	};

	~MatrixInterface()
	{
		clear_eigenstuff();
	}


//	void set_default_options();
//	void set_verbose_invert     (int yesno ){	verbose_invert      = yesno;	};
//	void set_verbose_diagonalize(int yesno ){	verbose_diagonalize = yesno;	};
//	void set_diagonalize_region (int region){	diagonalize_region  = region;	};
//	void set_diagonalize_numev  (int numev ){	diagonalize_numev   = numev;	};

	void clear();
	int size() { return A.rows(); };

	void MultDv(T *v, T *w)
	{
		//Count the number of times this method is called. Print every 100 times
		++Mv_counter;
		if(Mv_counter%100 == 0) { printf(" %li", Mv_counter); fflush(stdout); }

		Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1> > E_v(v,A.rows());
		Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1> > E_w(w,A.rows());
		E_w = A * E_v;
	};
	void cMultDv(T const *v, T *w) const
	{
		Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1> > E_v(v,A.rows());
		Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1> > E_w(w,A.rows());
		E_w = A * E_v;
	};

	T cMultvDw_real(T const *v, T const *w) const //THIS SHOULD BE COMBINED WITH cMultvDw!!!!!!!
	{
		T Dw[A.rows()];
		cMultDv(w, Dw);
		T result = T(0);
		for(int i=0; i<A.rows(); ++i)
			result += v[i]*Dw[i];
		return result;
	}

	T cMultvDw(T const *v, T const *w) const
	{
		T Dw[A.rows()];
		cMultDv(w, Dw);
		T result = T(0);
		for(int i=0; i<A.rows(); ++i)
			result += std::conj(v[i])*Dw[i];
		return result;
	}

	void MultOPv(T *v, T *w)
	{
		++OPv_counter;
		if(OPv_counter%100 == 0) { printf(" %li", OPv_counter); fflush(stdout); }
		Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1> > E_v(v,A.rows());
		Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1> > E_w(w,A.rows());
		E_w = OP_solver.solve(E_v);
	}

	void clear_eigenstuff()
	{
		if(Eigenvalues != NULL)
			delete[] Eigenvalues;
		if(Eigenvectors != NULL)
			delete[] Eigenvectors;
	};

	void Diagonalize();
//	void Diagonalize_shift_and_invert(T sig = (T)0.);
	void HumanReport(int num_print = 30) const;
//	void get_residuals(int nconv, T const *evals, T const *evecs, T *&residuals) const;
	void CalculateResidule(T *res) const;

	void LUdecomp();
	void solve_Mx_equals_b(T *x, T const *b);

//	void save_matrix_nc(char *filename);
//	void load_matrix_nc(char *filename);

};

template<typename T>
void MatrixInterface<T>::Diagonalize()
{
	clear_eigenstuff();

	std::cout << "Begin Diagonalization..." << std::endl;
	time_t start,end;
	time(&start);
	Mv_counter = 0;
	std::cout << "calls to MultDv:" << std::endl;
	fflush(stdout);

	if(num_request <= 0) num_request = A.rows()-1;
	else
		num_request=std::min(A.rows()-1, num_request);
	Eigenvalues = new T[num_request];
	if(compute_vecs)
		Eigenvectors = new T[num_request*A.rows()];

	ARWrapper<T> Wrapper;
	Wrapper.Solver.DefineParameters(A.rows(),num_request,this,&MatrixInterface<T>::MultDv,(char*)ARWrapper<T>::which.c_str());
	if(compute_vecs)
		num_converged = Wrapper.Solver.EigenValVectors(Eigenvectors,Eigenvalues);
	else
		num_converged = Wrapper.Solver.Eigenvalues(Eigenvalues);

	std::cout << std::endl;
	time(&end);
	diagonalization_time = difftime(end,start);
	num_Mv_calls = Mv_counter;

//	if(verbose_diagonalize)
//		eigenvalue_report();
};


/*
template<typename T>
void MatrixInterface<T>::Diagonalize_shift_and_invert(T sigma)
{
	clear_eigenstuff();

	std::cout << "Begin Diagonalization..." << std::endl;
	time_t start,end;
	time(&start);
	OPv_counter = 0;
	Mv_counter = 0;
	std::cout << "calls to MultOPv:" << std::endl;
	fflush(stdout);

	if(num_request <= 0) num_request = A.rows()-1;
	else
		num_request=std::min(A.rows()-1, num_request);
	Eigenvalues = new T[num_request];
	if(compute_vecs)
		Eigenvectors = new T[num_request*A.rows()];

//	EMatrix Identity(A.rows(), A.rows());
//	set_identity(Identity, A.rows());
//	OPinv = A - sigma*Identity;
//	OP_solver.compute(OPinv);
	ARWrapper<T,EMatrix> Wrapper;
//	Wrapper.Solver.DefineParameters(A.rows(),num_request,this,&MatrixInterface<T>::MultOPv,(char*)ARWrapper<T,EMatrix>::which2.c_str());
	Wrapper.Solver.DefineParameters(A.rows(),num_request,this,&MatrixInterface<T>::MultDv,(char*)ARWrapper<T,EMatrix>::which.c_str());
	Wrapper.Solver.ChangeShift(sigma);
	if(compute_vecs)
		num_converged = Wrapper.Solver.EigenValVectors(Eigenvectors,Eigenvalues);
	else
		num_converged = Wrapper.Solver.Eigenvalues(Eigenvalues);

	for(int i=0; i<num_converged; ++i)
		Eigenvalues[i] = ((T)1.)/Eigenvalues[i];

	std::cout << std::endl;
	time(&end);
	diagonalization_time = difftime(end,start);
	num_OPv_calls = OPv_counter;

//	if(verbose_diagonalize)
//		eigenvalue_report();
};
*/

template<typename T>
void MatrixInterface<T>::HumanReport(int num_print) const
{
	std::cout << "\nGenerating eigenvalue report..." << std::endl;

//	std::cout << "\tUsing ARPACK++ class ARSymStdEig" << std::endl;
//	std::cout << "\tReal symmetric eigenvalue problem: D*x - lambda*x" << std::endl;
	std::cout << "\tDimension of the system             : " << A.rows()				<< std::endl;
	std::cout << "\tNumber of 'requested' eigenvalues   : " << num_request			<< std::endl;
	std::cout << "\tNumber of 'converged' eigenvalues   : " << num_converged		<< std::endl;
	std::cout << "\tSeconds to diagonalize              : " << diagonalization_time	<< std::endl;
	std::cout << "\tNumber of calls to MultMv           : " << num_Mv_calls			<< std::endl << std::endl;

	num_print = std::min(num_print, num_converged);
	printf("     i:         eigenvalue\n");
	for(int i=0; i<num_print; ++i)
		printf("%6i:  % e \n", i, Eigenvalues[i]);
};
	
template<typename T>
void MatrixInterface<T>::CalculateResidule(T *res) const
{
	if(!compute_vecs)
	{
		printf("Cannot calculate the residules without the eigenvectors\n");
		return;
	}

	dbl rv, rw, iv, iw;
	const int Nvar = A.rows();
	std::complex<dbl> w[Nvar];

	for(int i=0; i<num_converged; ++i)
	{
		//  dbl *evec = modesAll[i].vec;
		cMultDv(&Eigenvectors[Nvar*i], w);
		std::complex<dbl> lambda(0., 0.);
		for(int j=0; j<Nvar; ++j)
		{   
			rv = std::real(Eigenvectors[Nvar*i+j]);
			iv = std::imag(Eigenvectors[Nvar*i+j]);
			rw = std::real(w[j]);
			iw = std::imag(w[j]);
			lambda += std::complex<dbl>(rv*rw+iv*iw, rv*iw-iv*rw);
		}
		res[i] = lambda - Eigenvalues[i];
	}
}

template<typename T>
void MatrixInterface<T>::LUdecomp()
{
	UMFsolver.compute(A);
	if(UMFsolver.info()!=Eigen::Success) 
	{
		printf("ERROR: Eigen::UmfPackLU decomposition failed\n");	
		exit(EXIT_FAILURE);
	}
}

template<typename T>
void MatrixInterface<T>::solve_Mx_equals_b(T *x, T const *b)
{
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VEC;
	Eigen::Map<      VEC> X(&x[0], A.rows());
	Eigen::Map<const VEC> B(&b[0], A.rows());
	X = UMFsolver.solve(B);
}


}

#endif //MATRIX_INTERFACE_H





