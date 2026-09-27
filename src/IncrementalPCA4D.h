// IncrementalPCA4D.hpp
//
// Incremental (online) Principal Component Analysis for exactly 4 input
// variables.
//
// - Data points are fed in one at a time via addSample().
// - Mean and covariance are updated incrementally (Welford-style running
//   covariance, numerically stable, O(1) memory, no need to store samples).
// - The PCA itself (eigenvectors = component coefficients / loadings,
//   eigenvalues = explained variance) is computed on demand from the
//   current covariance matrix using a classical Jacobi eigenvalue solver.
//
// Usage:
//   IncrementalPCA4D pca;
//   pca.addSample({1.0, 2.0, 0.5, 3.1});
//   pca.addSample({1.2, 1.9, 0.4, 3.0});
//   ...
//   auto components = pca.components();          // 4x4, columns = PCs
//   auto variance    = pca.explainedVariance();   // eigenvalues, descending
//   auto ratio       = pca.explainedVarianceRatio();
//
// Single header, no dependencies beyond <array>, <cmath>, <algorithm>.

#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

template<int DIM>
class IncrementalPCA4D {
public:
   //static constexpr int N = 4; // fixed number of variables
    static constexpr int N = DIM; // fixed number of variables

    using Vector4 = std::array<double, N>;
    using Matrix4 = std::array<std::array<double, N>, N>;

    IncrementalPCA4D() {
        mean_.fill(0.0);
        for (auto& row : M2_) row.fill(0.0);
    }

    /// Reset the accumulator to the empty state.
    void reset() {
        n_ = 0;
        mean_.fill(0.0);
        for (auto& row : M2_) row.fill(0.0);
    }

    /// Add one new 4-dimensional sample. O(N^2) = O(16) per call.
    void addSample(const Vector4& x) {
        ++n_;
        Vector4 deltaOld{};
        for (int k = 0; k < N; ++k) deltaOld[k] = x[k] - mean_[k];
        for (int k = 0; k < N; ++k) mean_[k] += deltaOld[k] / static_cast<double>(n_);
        Vector4 deltaNew{};
        for (int k = 0; k < N; ++k) deltaNew[k] = x[k] - mean_[k];
        // Running sum of outer products of deviations (Welford, multivariate form).
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                M2_[i][j] += deltaOld[i] * deltaNew[j];
    }

    /// Number of samples seen so far.
    long long sampleCount() const { return n_; }

    /// Current running mean.
    Vector4 mean() const { return mean_; }

    /// Current sample covariance matrix (unbiased, divides by n-1).
    /// Returns a zero matrix if fewer than 2 samples have been added.
    Matrix4 covariance() const {
        Matrix4 cov{};
        if (n_ < 2) {
            for (auto& row : cov) row.fill(0.0);
            return cov;
        }
        const double denom = static_cast<double>(n_ - 1);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                cov[i][j] = M2_[i][j] / denom;
        return cov;
    }

    Matrix4 correlation() const {
        Matrix4 cov{};
        if (n_ < 2) {
            for (auto& row : cov) row.fill(0.0);
            return cov;
        }
        const double denom = static_cast<double>(n_ - 1);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                cov[i][j] = M2_[i][j] / denom;

        Matrix4 cor{};
        for (int i = 0; i < N; ++i)
	{
            for (int j = 0; j < N; ++j)
		if(i != j)
                    cor[i][j] = cov[i][j]/std::sqrt(cov[i][i]*cov[j][j]);
		else
	            cor[i][i] = 1;
	}
        return cor;
    }

    /// Principal component coefficients (loadings).
    /// Returned matrix: components()[i][k] is the coefficient of original
    /// variable i in principal component k. Columns are unit-length
    /// eigenvectors of the covariance matrix, sorted by descending
    /// explained variance (component 0 = largest variance).
    Matrix4 components() const {
        Matrix4 V{};
        Vector4 lambda{};
        computePCA(V, lambda);
        return V;
    }

    /// Eigenvalues of the covariance matrix (= variance explained by each
    /// principal component), sorted descending, matching components().
    Vector4 explainedVariance() const {
        Matrix4 V{};
        Vector4 lambda{};
        computePCA(V, lambda);
        return lambda;
    }

    /// Same as explainedVariance(), normalized to fractions summing to 1.
    Vector4 explainedVarianceRatio() const {
        Vector4 lambda = explainedVariance();
        double total = 0.0;
        for (double v : lambda) total += v;
        Vector4 ratio{};
        for (int i = 0; i < N; ++i) ratio[i] = (total > 0.0) ? lambda[i] / total : 0.0;
        return ratio;
    }

    /// Project a (mean-centered) sample onto the principal components,
    /// i.e. compute its PCA scores.
    Vector4 transform(const Vector4& x) const {
        Matrix4 V = components();
        Vector4 centered{};
        for (int i = 0; i < N; ++i) centered[i] = x[i] - mean_[i];
        Vector4 scores{};
        for (int k = 0; k < N; ++k) {
            double s = 0.0;
            for (int i = 0; i < N; ++i) s += centered[i] * V[i][k];
            scores[k] = s;
        }
        return scores;
    }

private:
    long long n_ = 0;
    Vector4 mean_{};
    Matrix4 M2_{}; // running sum of outer products of deviations

    /// Runs a Jacobi eigenvalue decomposition of the current covariance
    /// matrix and returns eigenvectors (columns of V) / eigenvalues,
    /// sorted by descending eigenvalue.
    void computePCA(Matrix4& V, Vector4& lambda) const {
        Matrix4 A = covariance();
        jacobiEigenDecomposition(A, V, lambda);
        sortDescending(V, lambda);
    }

    /// Classical (cyclic) Jacobi eigenvalue algorithm for a real symmetric
    /// NxN matrix. On return, V's columns are orthonormal eigenvectors and
    /// lambda holds the corresponding eigenvalues (unsorted).
    static void jacobiEigenDecomposition(Matrix4 A, Matrix4& V, Vector4& lambda,
                                          int maxSweeps = 100, double tol = 1e-14) {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                V[i][j] = (i == j) ? 1.0 : 0.0;

        for (int sweep = 0; sweep < maxSweeps; ++sweep) {
            double off = 0.0;
            for (int i = 0; i < N; ++i)
                for (int j = i + 1; j < N; ++j)
                    off += A[i][j] * A[i][j];
            if (off < tol) break;

            for (int p = 0; p < N - 1; ++p) {
                for (int q = p + 1; q < N; ++q) {
                    if (std::abs(A[p][q]) < 1e-300) continue;

                    double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                    double t = (theta >= 0.0 ? 1.0 : -1.0) /
                               (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                    double c = 1.0 / std::sqrt(t * t + 1.0);
                    double s = t * c;

                    double app = A[p][p], aqq = A[q][q], apq = A[p][q];
                    A[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                    A[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                    A[p][q] = 0.0;
                    A[q][p] = 0.0;

                    for (int i = 0; i < N; ++i) {
                        if (i == p || i == q) continue;
                        double aip = A[i][p], aiq = A[i][q];
                        A[i][p] = c * aip - s * aiq;
                        A[p][i] = A[i][p];
                        A[i][q] = s * aip + c * aiq;
                        A[q][i] = A[i][q];
                    }
                    for (int i = 0; i < N; ++i) {
                        double vip = V[i][p], viq = V[i][q];
                        V[i][p] = c * vip - s * viq;
                        V[i][q] = s * vip + c * viq;
                    }
                }
            }
        }
        for (int i = 0; i < N; ++i) lambda[i] = A[i][i];
    }

    /// Sort eigenpairs (columns of V, entries of lambda) by descending
    /// eigenvalue, in place.
    static void sortDescending(Matrix4& V, Vector4& lambda) {
        std::array<int, N> idx{};
	for(int i = 0; i < N; i++) idx[i] = i;
        std::sort(idx.begin(), idx.end(),
                  [&](int a, int b) { return lambda[a] > lambda[b]; });

        Matrix4 Vsorted{};
        Vector4 lsorted{};
        for (int k = 0; k < N; ++k) {
            lsorted[k] = lambda[idx[k]];
            for (int i = 0; i < N; ++i) Vsorted[i][k] = V[i][idx[k]];
        }
        V = Vsorted;
        lambda = lsorted;
    }
};
