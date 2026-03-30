#ifndef PCA_H
#define PCA_H

#include <vector>

namespace PCA {
    extern int N;

    std::vector<double> pca(const std::vector<double>& X, const std::vector<std::vector<double>>& eVec);
    void updateLgs(const std::vector<double>& X, double Y, int offset = 0);
    void corr(const std::vector<double>& X, int offset = 0);
    void solveLgs(int offset = 0);
}

#endif
