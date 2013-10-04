#ifndef CCS_MATRIX_H
#define CCS_MATRIX_H

#include "../Resources/std_include.h"
#include <Eigen/Sparse>

template <typename CCSint=int>
class CCSmatrix
{
	typedef Eigen::Matrix<CCSint, Eigen::Dynamic, 1> intvec;
	typedef Eigen::Matrix<dbl,    Eigen::Dynamic, 1> dblvec;
	
public:
	CCSint n_rows;
	CCSint n_cols;
	intvec Ap;
	intvec Ai;
	dblvec Ax;
	
	CCSmatrix();
	CCSmatrix(Eigen::SparseMatrix<dbl> const &A);
	
	void SetMatrix(Eigen::SparseMatrix<dbl> const &A);
	void Av(dbl const *v, dbl *w) const;
	void Av(Eigen::VectorXd const &V, Eigen::VectorXd &W) const;
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename CCSint>
CCSmatrix<CCSint>::CCSmatrix()
{}

template <typename CCSint>
CCSmatrix<CCSint>::CCSmatrix(Eigen::SparseMatrix<dbl> const &A)
{
	SetMatrix(A);
}

template <typename CCSint>
void CCSmatrix<CCSint>::SetMatrix(Eigen::SparseMatrix<dbl> const &A)
{
	n_rows = A.rows();
	n_cols = A.cols();

	int innerSize = A.innerSize();
	int outerSize = A.outerSize();
	int NnonZeros = A.nonZeros();

	Ap.resize(n_cols+1);
	Ai.resize(NnonZeros);
	Ax.resize(NnonZeros);

	for(CCSint ii = 0; ii<n_cols+1; ++ii)
		Ap[ii] = (CCSint)A.outerIndexPtr()[ii];

	for(CCSint ii = 0; ii<NnonZeros; ++ii)
	{
		Ai[ii] = (CCSint)A.innerIndexPtr()[ii];
		Ax[ii] = A.valuePtr()[ii];
	}
}

//Calculate w=Av
template <typename CCSint>
void CCSmatrix<CCSint>::Av(dbl const *v, dbl *w) const
{
	int r, c, i;
	for(r=0; r<n_rows; ++r)
		w[r]=0.;
	for(c=0; c<n_cols; ++c)
		for(i=Ap[c]; i<Ap[c+1]; ++i)
			w[Ai[i]] += v[c]*Ax[i];
}

template <typename CCSint>
void CCSmatrix<CCSint>::Av(Eigen::VectorXd const &V, Eigen::VectorXd &W) const
{
	assert(V.size() == n_cols);
	assert(W.size() == n_rows);
	Av(V.data(), W.data());
}



#endif //CCS_MATRIX_H
