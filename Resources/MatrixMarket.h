#ifndef MATRIX_MARKET_H
#define MATRIX_MARKET_H


#include <Eigen/Core>
#include <Eigen/Sparse>


class MatrixMarket
{
	static const string comment;

public:
	template <typename Derived>
	static void WriteDenseDbl	(char *filename, Eigen::MatrixBase<Derived> &A);
	static void WriteSparseDbl	(char *filename, Eigen::SparseMatrix<dbl> &A, bool is_symmetric=false);
	static void ReadVectorDbl	(char *filename, Eigen::VectorXd &V);
};



template <typename Derived>
void MatrixMarket::WriteDenseDbl(char *filename, Eigen::MatrixBase<Derived> &A)
{
	bool is_symmetric = false;
	FILE *outfile = fopen(filename, "w");

	fprintf(outfile, "%%%%MatrixMarket matrix array real %s\n", is_symmetric?"symmetric":"general");
	int nrows, ncols;
	nrows = A.rows();
	ncols = A.cols();
	fprintf(outfile, "%i %i\n", nrows, ncols);

	for(int c=0; c<ncols; ++c)
		for(int r=0; r<nrows; ++r)
			fprintf(outfile, "% 18.16e\n", A(r,c));
}



void MatrixMarket::WriteSparseDbl(char *filename, Eigen::SparseMatrix<dbl> &A, bool is_symmetric)
{
	FILE *outfile = fopen(filename, "w");

	fprintf(outfile, "%%%%MatrixMarket matrix coordinate real %s\n", is_symmetric?"symmetric":"general");
	fprintf(outfile, "%s", comment.c_str());
	int nrows, ncols, nnonzeros;
	nrows = A.rows();
	ncols = A.cols();
	if(is_symmetric)
	{
		assert(nrows == ncols);
		nnonzeros = 0;
		for (int k=0; k<A.outerSize(); ++k)
			for (Eigen::SparseMatrix<dbl>::InnerIterator it(A,k); it; ++it)
				if(it.row() >= it.col())
					++nnonzeros;
	}else{
		nnonzeros = A.nonZeros();
	}
	fprintf(outfile, "%i %i %i\n", nrows, ncols, nnonzeros);

	int counter = 0;
	for (int k=0; k<A.outerSize(); ++k)
		for (Eigen::SparseMatrix<dbl>::InnerIterator it(A,k); it; ++it)
		{
			if(!is_symmetric || it.row() >= it.col())
			{
				++counter;
				fprintf(outfile, "%10i %10i % 18.16e\n", it.row()+1, it.col()+1, it.value()); //The +1 is needed because MatrixMarket indexing starts at 1 not 0
			}
		}

	assert(counter == nnonzeros);

	fflush(outfile);
	fclose(outfile);

}

void MatrixMarket::ReadVectorDbl(char *filename, Eigen::VectorXd &V)
{
	string line;
	ifstream myfile (filename);

	if (myfile.is_open())
	{
		getline(myfile, line);
		//Check first line
		//...

		//Skip the comments and get the important first line
		getline(myfile, line);
		while (line[0]=='%')
			getline(myfile, line);

		//Read the # of rows and columns
		int nrows, ncols;
		istringstream buffer(line);
		buffer >> nrows >> ncols;
		
		V = Eigen::VectorXd::Zero(nrows);

		//Read the data
		int index = 0;
		while ( getline (myfile,line) )
		{
			V[index] = atof(line.c_str());
			++index;
		}

		assert(index == V.size());
		myfile.close();
	}
	else cout << "Unable to open file\n"; 
}


const string MatrixMarket::comment = "\
%=================================================================================\n\
%\n\
% This ASCII file represents a sparse MxN matrix with L \n\
% nonzeros in the following Matrix Market format:\n\
%\n\
% +----------------------------------------------+\n\
% |%%MatrixMarket matrix coordinate real general | <--- header line\n\
% |%                                             | <--+\n\
% |% comments                                    |    |-- 0 or more comment lines\n\
% |%                                             | <--+         \n\
% |    M  N  L                                   | <--- rows, columns, entries\n\
% |    I1  J1  A(I1, J1)                         | <--+\n\
% |    I2  J2  A(I2, J2)                         |    |\n\
% |    I3  J3  A(I3, J3)                         |    |-- L lines\n\
% |        . . .                                 |    |\n\
% |    IL JL  A(IL, JL)                          | <--+\n\
% +----------------------------------------------+   \n\
%\n\
% Indices are 1-based, i.e. A(1,1) is the first element.\n\
%\n\
%=================================================================================\n";




#endif //MATRIX_MARKET_H














