// incremental_logistic_regression.cpp
//
// Model:   logit(P(T=1)) = b0 + b1*D + b2*C + b3*(D*C)
// where T is a binary label, C in {0,1}, D is a numeric feature.
//
// b3 is exactly the difference between the log-odds slope of D for
// C=1 and the log-odds slope of D for C=0. Its Wald z-statistic /
// p-value tells you whether the two "A" slopes differ significantly
// (the logistic-regression analogue of the linear-regression case).
//
// IMPORTANT DIFFERENCE VS. LINEAR REGRESSION:
// Ordinary least squares has a closed form, so we could keep only a
// few running sums (X'X, X'y) and fit in O(1) per new point. Logistic
// regression has NO closed form: the maximum-likelihood fit requires
// an iterative solver (IRLS / Newton-Raphson) whose weights depend
// non-linearly on the *current* coefficient estimate. That means we
// cannot compress the data into fixed-size sufficient statistics.
//
// So here "incremental" means: adding a point is O(1) (just store
// it), and fit() re-runs IRLS over all stored points whenever you
// ask for the current estimate. This is the standard, correct way to
// do it -- it is just not O(1) per fit the way linear regression was.
// (A fully O(1)-per-point online approximation exists -- online
// Newton step / online IRLS with a running Hessian approximation --
// but it trades exactness for speed and is not what most people mean
// when they ask for "the" logistic regression fit.)
//
// Compile:  g++ -O2 -std=c++17 incremental_logistic_regression.cpp -o logregtest
// Run:      ./logregtest

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

// ---------- tiny fixed-size (4x4) linear algebra ----------

using Mat4 = std::array<std::array<double, 4>, 4>;
using Vec4 = std::array<double, 4>;

static Mat4 invert4(Mat4 a) {
    Mat4 inv{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            inv[i][j] = (i == j) ? 1.0 : 0.0;

    for (int col = 0; col < 4; ++col) {
        int pivotRow = col;
        double best = std::fabs(a[col][col]);
        for (int r = col + 1; r < 4; ++r) {
            if (std::fabs(a[r][col]) > best) {
                best = std::fabs(a[r][col]);
                pivotRow = r;
            }
        }
        if (best < 1e-12)
	{
            //throw std::runtime_error("Matrix is singular (not enough / degenerate data yet).");
	    std::cerr << "Abort: Matrix is singular (not enough / degenerate data yet)." << std::endl;
	    std::exit(-2);
	}

        std::swap(a[col], a[pivotRow]);
        std::swap(inv[col], inv[pivotRow]);

        double piv = a[col][col];
        for (int j = 0; j < 4; ++j) { a[col][j] /= piv; inv[col][j] /= piv; }

        for (int r = 0; r < 4; ++r) {
            if (r == col) continue;
            double factor = a[r][col];
            if (factor == 0.0) continue;
            for (int j = 0; j < 4; ++j) {
                a[r][j] -= factor * a[col][j];
                inv[r][j] -= factor * inv[col][j];
            }
        }
    }
    return inv;
}

// standard normal CDF, for Wald z-test p-values
static double normalCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

static double sigmoid(double eta) {
    if (eta >= 0) {
        double e = std::exp(-eta);
        return 1.0 / (1.0 + e);
    } else {
        double e = std::exp(eta);
        return e / (1.0 + e);
    }
}

// ---------- incremental (point-at-a-time) logistic regression ----------

struct Observation {
    double D, C, T; // T must be 0.0 or 1.0
};

struct FitResult {
    Vec4 beta{};   // b0, b1(D), b2(C), b3(D*C)  <-- b3 is the slope difference
    Vec4 se{};     // standard errors (from inverse Fisher information)
    Vec4 zstat{};
    Vec4 pvalue{}; // two-sided Wald test p-value
    double logLikelihood = 0.0;
    int iterations = 0;
    bool converged = false;
};

class IncrementalLogisticRegression {
public:
    // Add one observation. C should be 0.0 or 1.0, T should be 0.0 or 1.0.
    void addPoint(double D, double C, double T) {
        data_.push_back({D, C, T});
    }

    long count() const { return static_cast<long>(data_.size()); }

    // Re-fits by maximum likelihood (IRLS / Newton-Raphson) over all
    // points added so far.
    FitResult fit(bool debug = false, int maxIter = 100, double tol = 1e-9) const {
        long n = count();
        if (n < 5)
	{
            //throw std::runtime_error("Need at least 5 points to attempt a fit.");
	    std::cerr << "Need at least 5 points to attempt a fit." << std::endl;
	    std::exit(-1);
	}

        Vec4 beta{0.0, 0.0, 0.0, 0.0};
        const double eps = 1e-9; // clamp to avoid weights collapsing to 0

        FitResult r;
        Mat4 XtWXinv{};

        for (int iter = 0; iter < maxIter; ++iter) {
            Mat4 XtWX{};
            Vec4 XtWz{};

            for (const auto& obs : data_) {
                Vec4 x = {1.0, obs.D, obs.C, obs.D * obs.C};
                double eta = 0.0;
                for (int k = 0; k < 4; ++k) eta += beta[k] * x[k];
                double p = sigmoid(eta);
                p = std::min(std::max(p, eps), 1.0 - eps);
                double w = p * (1.0 - p);
                double z = eta + (obs.T - p) / w; // working response

                for (int i = 0; i < 4; ++i) {
                    XtWz[i] += x[i] * w * z;
                    for (int j = 0; j < 4; ++j)
                        XtWX[i][j] += x[i] * w * x[j];
                }
            }

            XtWXinv = invert4(XtWX);
            Vec4 newBeta{};
            for (int i = 0; i < 4; ++i) {
                double s = 0.0;
                for (int j = 0; j < 4; ++j) s += XtWXinv[i][j] * XtWz[j];
                newBeta[i] = s;
            }

            double maxDelta = 0.0;
            for (int i = 0; i < 4; ++i)
                maxDelta = std::max(maxDelta, std::fabs(newBeta[i] - beta[i]));

            beta = newBeta;
            r.iterations = iter + 1;

	    if(debug)
	    {
		    std::cerr << "Finished iteration " << iter+1 << " error=" << maxDelta << std::endl;
		    std::cerr << "=> beta:";
		    for (int i = 0; i < 4; ++i)
			    std::cerr << " " << beta[i];
		    std::cerr << std::endl;
	    }

            if (maxDelta < tol) { r.converged = true; break; }
        }

        // Log-likelihood and Wald standard errors at the converged beta.
        double ll = 0.0;
        for (const auto& obs : data_) {
            Vec4 x = {1.0, obs.D, obs.C, obs.D * obs.C};
            double eta = 0.0;
            for (int k = 0; k < 4; ++k) eta += beta[k] * x[k];
            double p = sigmoid(eta);
            p = std::min(std::max(p, eps), 1.0 - eps);
            ll += obs.T * std::log(p) + (1.0 - obs.T) * std::log(1.0 - p);
        }

        r.beta = beta;
        r.logLikelihood = ll;
        for (int i = 0; i < 4; ++i) {
            r.se[i] = std::sqrt(XtWXinv[i][i]); // inverse observed Fisher information
            r.zstat[i] = (r.se[i] > 0.0) ? r.beta[i] / r.se[i] : 0.0;
            r.pvalue[i] = 2.0 * (1.0 - normalCDF(std::fabs(r.zstat[i])));
        }
        return r;
    }

private:
    std::vector<Observation> data_;
};

// ---------- demo / usage ----------

static void printFit(const FitResult& r, long n, std::ostream& out = std::cerr) {
    static const char* names[4] = {"Intercept (b0)", "D slope (b1)",
                                    "C offset (b2)", "D*C interaction (b3)"};
    out << "n = " << n
              << ", converged = " << (r.converged ? "yes" : "no")
              << " (" << r.iterations << " iterations)"
              << ", logLik = " << r.logLikelihood << "\n";
    out.setf(std::ios::fixed);
    out.precision(5);
    for (int i = 0; i < 4; ++i) {
        out << "  " << names[i]
                  << ": beta=" << r.beta[i]
                  << "  se=" << r.se[i]
                  << "  z=" << r.zstat[i]
                  << "  p=" << r.pvalue[i]
                  << "\n";
    }
    out << "  -> Log-odds slope for C=false : " << r.beta[1] << "\n";
    out << "  -> Log-odds slope for C=true  : " << (r.beta[1] + r.beta[3]) << "\n";
    out << "  -> Slope difference (A_true - A_false) = b3 = " << r.beta[3]
              << ", p = " << r.pvalue[3];
    out << (r.pvalue[3] < 0.05 ? "  => significant at alpha=0.05\n"
                                      : "  => not significant at alpha=0.05\n");
}

/*
int main() {
    IncrementalLogisticRegression reg;

    // ---- Example: synthetic data, added ONE POINT AT A TIME ----
    // True model, C=false: logit(P) = -3.0 + 1.0*D
    // True model, C=true : logit(P) = -3.0 + 3.0*D   (different slope!)
    // Replace this generator with your real incremental data source --
    // just call reg.addPoint(D, C, T) whenever a new observation arrives.

    std::mt19937 rng(42); // fixed seed for reproducible demo output
    std::uniform_real_distribution<double> dDist(0.2, 3.0);
    std::uniform_real_distribution<double> uDist(0.0, 1.0);

    const double b0 = -3.0, bFalse = 1.0, bTrue = 3.0;
    const int perGroup = 60;

    std::vector<Observation> pool;
    for (int i = 0; i < perGroup; ++i) {
        double D = dDist(rng);
        double etaF = b0 + bFalse * D;
        double pF = sigmoid(etaF);
        double T_F = (uDist(rng) < pF) ? 1.0 : 0.0;
        pool.push_back({D, 0.0, T_F});

        D = dDist(rng);
        double etaT = b0 + bTrue * D;
        double pT = sigmoid(etaT);
        double T_T = (uDist(rng) < pT) ? 1.0 : 0.0;
        pool.push_back({D, 1.0, T_T});
    }
    // shuffle so points arrive in a realistic mixed order
    std::shuffle(pool.begin(), pool.end(), rng);

    for (const auto& obs : pool) {
        reg.addPoint(obs.D, obs.C, obs.T);

        if (reg.count() % 20 == 0) {
            std::cout << "---- after " << reg.count() << " points ----\n";
            try {
                printFit(reg.fit(), reg.count());
            } catch (const std::exception& e) {
                std::cout << "  (" << e.what() << ")\n";
            }
        }
    }

    std::cout << "\n==== FINAL RESULT ====\n";
    printFit(reg.fit(), reg.count());

    return 0;
}
*/
