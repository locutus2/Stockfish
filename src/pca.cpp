#include <iostream>

#include "pca.h"
#include "misc.h"

namespace PCA {

int N = 0;

// https://www.emathhelp.net/calculators/linear-algebra/eigenvalue-and-eigenvector-calculator
std::vector<double> pca(const std::vector<double>& X, const std::vector<std::vector<double>>& eVec)
{
        const int n = int(X.size());
        std::vector<double> v;
        for(int i = 0; i < int(eVec.size()); i++)
        {
            double s = 0;
            for(int j = 0; j < n; j++)
                    s += X[j] * eVec[i][j];
            v.push_back(s);
        }
        return v;
}

void updateLgs(const std::vector<double>& X, double Y, int offset)
{
        const int n = int(X.size());
        N = n;
        for(int i = 0; i < n; i++)
           for(int j = 0; j < n; j++)
           {
		   Stockfish::dbg_sum_of(X[i] * X[j], offset+10*i+j);
           }

        for(int i = 0; i < n; i++)
        {
		Stockfish::dbg_sum_of(X[i] * Y, offset+100+i);
        }
}

void corr(const std::vector<double>& X, int offset)
{
        const int n = int(X.size());
        for(int i = 0; i < n; i++)
           for(int j = 0; j < n; j++)
           {
		   Stockfish::dbg_correl_of(X[i], X[j], offset+10*i+j);
           }
}

void solveLgs(int offset)
{
        int n = N;

	std::vector<std::vector<long double>> A(n, std::vector<long double>(n, 0));
	std::vector<std::vector<long double>> Ainv(n, std::vector<long double>(n, 0));

	for(int i = 0; i < n; i++)
		Ainv[i][i] = 1;

	for(int i = 0; i < n; i++)
	     for(int j = 0; j < n; j++)
		     A[i][j] = Stockfish::dbg_get_sum_of(offset+10*i+j);

	// forward
	for(int i = 0; i < n; i++)
	{
		assert(A[i][i] != 0);
	     for(int k = i+1; k < n; k++)
	     {
	          for(int j = i+1; j < n; j++)
		     A[k][j] -= A[i][j] * A[k][i] / A[i][i];
	          for(int j = 0; j < n; j++)
		     Ainv[k][j] -= Ainv[i][j] * A[k][i] / A[i][i];

		  A[k][i] = 0;
	     }
	}

	// backward
	for(int i = n-1; i >= 0; i--)
	{
		assert(A[i][i] != 0);
		
	        for(int j = i+1; j < n; j++)
			A[i][j] /= A[i][i];
	        for(int j = 0; j < n; j++)
			Ainv[i][j] /= A[i][i];

		A[i][i] = 1;
	        for(int k = i-1; k >= 0; k--)
	        { 
	            for(int j = 0; j < n; j++)
			    Ainv[k][j] -= Ainv[i][j] * A[k][i];
		    A[k][i] = 0;
	        }
	}

        std::vector<long double> beta(n, 0);
        for(int i = 0; i < n; i++)
        {
             for(int j = 0; j < n; j++)
		     beta[i] += Ainv[i][j] * Stockfish::dbg_get_sum_of(offset + 100 + j);
        }

        std::cerr << "I:" << std::endl;
        for(int i = 0; i < n; i++)
        {
             for(int j = 0; j < n; j++)
		     std::cerr << A[i][j] << "\t";
             std::cerr << std::endl;
        }
        std::cerr << "Inverse XtX:" << std::endl;
        for(int i = 0; i < n; i++)
        {
             for(int j = 0; j < n; j++)
		     std::cerr << Ainv[i][j] << "\t";
             std::cerr << std::endl;
        }
        std::cerr << "XtY:" << std::endl;
        for(int i = 0; i < n; i++)
        {
	    std::cerr << Stockfish::dbg_get_sum_of(offset + 100 + i) << "\t";
        }
        std::cerr << std::endl;

        std::cerr << "beta:" << std::endl;
        for(int i = 0; i < n; i++)
        {
		     std::cerr << beta[i] << "\t";
        }
        std::cerr << std::endl;
}

}

