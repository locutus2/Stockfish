// incremental_logistic_regression_quad.cpp
//
// Model (extended with quadratic terms):
//
//   logit(P(T=1)) = b0 + b1*D + b2*C + b3*(D*C) + b4*(D*D) + b5*(D*D*C)
//
// where T is a binary label, C in {0,1}, D is a numeric feature.
//
// This lets the log-odds curve be a full quadratic in D, with BOTH the
// linear slope and the curvature allowed to differ between the C=0
// and C=1 groups:
//   - b3 (D*C)   is the difference in the *linear* slope between groups.
//   - b5 (D*D*C) is the difference in the *curvature* (quadratic term)
//                between groups.
// Each still gets its own Wald z-test / p-value, same as before.
//
// Since b3 and b5 together describe "does the whole curve differ
// between groups", a joint test is often more meaningful than looking
// at each Wald test individually (they can each look non-significant
// while jointly being significant, or vice versa). For that this file
// also fits the nested REDUCED model
//
//   logit(P(T=1)) = b0 + b1*D + b2*C + b4*(D*D)
//
// (same shape curve for both groups, only a different intercept/
// baseline offset) and reports a likelihood-ratio test comparing it
// to the full model above. Dropping exactly 2 parameters (b3, b5)
// means the LR statistic follows a chi-square distribution with
// df = 2, whose CDF has the closed form P(X <= x) = 1 - exp(-x/2),
// so the LRT p-value can be computed exactly without any special
// functions: p = exp(-LR/2).
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
//
// Compile:  g++ -O2 -std=c++17 incremental_logistic_regression_quad.cpp -o logregtest
// Run:      ./logregtest

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

// ---------- tiny generic (NxN) linear algebra ----------
// Templated on size so we can invert both the 6x6 (full model) and
// the 4x4 (reduced model) matrices with one implementation.

template <std::size_t N>
using Mat = std::array<std::array<double, N>, N>;
template <std::size_t N>
using Vec = std::array<double, N>;

template <std::size_t N>
static Mat<N> invert(Mat<N> a) {
    Mat<N> inv{};
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            inv[i][j] = (i == j) ? 1.0 : 0.0;

    for (std::size_t col = 0; col < N; ++col) {
        std::size_t pivotRow = col;
        double best = std::fabs(a[col][col]);
        for (std::size_t r = col + 1; r < N; ++r) {
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
        for (std::size_t j = 0; j < N; ++j) { a[col][j] /= piv; inv[col][j] /= piv; }

        for (std::size_t r = 0; r < N; ++r) {
            if (r == col) continue;
            double factor = a[r][col];
            if (factor == 0.0) continue;
            for (std::size_t j = 0; j < N; ++j) {
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
    int8_t D, C, T; // T must be 0 or 1; upgraded from int8_t so D*D / D*D*C don't overflow the *value range* headroom when D is a few dozen or more
};

template <std::size_t N>
struct FitResultT {
    Vec<N> beta{};   // fitted coefficients
    Vec<N> se{};     // standard errors (from inverse Fisher information)
    Vec<N> zstat{};
    Vec<N> pvalue{}; // two-sided Wald test p-value
    double logLikelihood = 0.0;
    int iterations = 0;
    bool converged = false;
};

// Full model:    [1, D, C, D*C, D*D, D*D*C]   (6 parameters)
using FitResult = FitResultT<6>;
// Reduced model: [1, D, C, D*D]               (4 parameters, no interaction terms)
using FitResultReduced = FitResultT<4>;

struct LikelihoodRatioTest {
    double llFull = 0.0;
    double llReduced = 0.0;
    double lrStatistic = 0.0; // 2 * (llFull - llReduced)
    int df = 2;               // b3 and b5 dropped
    double pvalue = 1.0;      // exact, since df == 2: p = exp(-LR/2)
};

class IncrementalLogisticRegression {
public:
    // Add one observation. C should be 0 or 1, T should be 0 or 1.
    void addPoint(int8_t D, int8_t C, int8_t T) {
        data_.push_back({D, C, T});
    }

    long count() const { return static_cast<long>(data_.size()); }

    // Full model: logit(P) = b0 + b1*D + b2*C + b3*(D*C) + b4*(D*D) + b5*(D*D*C)
    FitResult fit(bool debug = false, int maxIter = 100, double tol = 1e-8) const {
        return fitGeneric<6>(
            [](const Observation& obs) -> Vec<6> {
                double D = obs.D, C = obs.C;
                return {1.0, D, C, D * C, D * D, D * D * C};
            },
            "full", debug, maxIter, tol);
    }

    // Reduced (nested) model: logit(P) = b0 + b1*D + b2*C + b4*(D*D)
    // Same curve shape for both groups; only the intercept differs by C.
    FitResultReduced fitReduced(bool debug = false, int maxIter = 100, double tol = 1e-8) const {
        return fitGeneric<4>(
            [](const Observation& obs) -> Vec<4> {
                double D = obs.D, C = obs.C;
                return {1.0, D, C, D * D};
            },
            "reduced", debug, maxIter, tol);
    }

    // Likelihood-ratio test of "do the two group curves differ at all"
    // (H0: b3 = 0 and b5 = 0, i.e. the reduced model holds).
    LikelihoodRatioTest likelihoodRatioTest(bool debug = false) const {
        FitResult full = fit(debug);
        FitResultReduced reduced = fitReduced(debug);

        LikelihoodRatioTest t;
        t.llFull = full.logLikelihood;
        t.llReduced = reduced.logLikelihood;
        t.lrStatistic = 2.0 * (t.llFull - t.llReduced);
        if (t.lrStatistic < 0.0) t.lrStatistic = 0.0; // guard against tiny numerical noise
        t.df = 2;
        t.pvalue = std::exp(-t.lrStatistic / 2.0); // exact chi-square(df=2) tail
        return t;
    }

private:
    // Shared IRLS core. FeatureFn maps an Observation to its design
    // vector (Vec<N>) for whichever model (full or reduced) is being fit.
    template <std::size_t N, typename FeatureFn>
    FitResultT<N> fitGeneric(FeatureFn features, const char* label, bool debug,
                              int maxIter, double tol) const {
        long n = count();
        if (n < 5)
	{
            //throw std::runtime_error("Need at least 5 points to attempt a fit.");
	    std::cerr << "Need at least 5 points to attempt a fit." << std::endl;
	    std::exit(-1);
	}

        Vec<N> beta{};
        const double eps = 1e-9; // clamp to avoid weights collapsing to 0

        FitResultT<N> r;
        Mat<N> XtWXinv{};

        for (int iter = 0; iter < maxIter; ++iter) {
            Mat<N> XtWX{};
            Vec<N> XtWz{};

            for (const auto& obs : data_) {
                Vec<N> x = features(obs);
                double eta = 0.0;
                for (std::size_t k = 0; k < N; ++k) eta += beta[k] * x[k];
                double p = sigmoid(eta);
                p = std::min(std::max(p, eps), 1.0 - eps);
                double w = p * (1.0 - p);
                double z = eta + (obs.T - p) / w; // working response

                for (std::size_t i = 0; i < N; ++i) {
                    XtWz[i] += x[i] * w * z;
                    for (std::size_t j = 0; j < N; ++j)
                        XtWX[i][j] += x[i] * w * x[j];
                }
            }

            XtWXinv = invert<N>(XtWX);
            Vec<N> newBeta{};
            for (std::size_t i = 0; i < N; ++i) {
                double s = 0.0;
                for (std::size_t j = 0; j < N; ++j) s += XtWXinv[i][j] * XtWz[j];
                newBeta[i] = s;
            }

            double maxDelta = 0.0;
            for (std::size_t i = 0; i < N; ++i)
                maxDelta = std::max(maxDelta, std::fabs(newBeta[i] - beta[i]));

            beta = newBeta;
            r.iterations = iter + 1;

	    if(debug)
	    {
		    std::cerr << "[" << label << "] Finished iteration " << iter+1 << " error=" << maxDelta << std::endl;
		    std::cerr << "=> beta:";
		    for (std::size_t i = 0; i < N; ++i)
			    std::cerr << " " << beta[i];
		    std::cerr << std::endl;
	    }

            if (maxDelta < tol) { r.converged = true; break; }
        }

        // Log-likelihood and Wald standard errors at the converged beta.
        double ll = 0.0;
        for (const auto& obs : data_) {
            Vec<N> x = features(obs);
            double eta = 0.0;
            for (std::size_t k = 0; k < N; ++k) eta += beta[k] * x[k];
            double p = sigmoid(eta);
            p = std::min(std::max(p, eps), 1.0 - eps);
            ll += obs.T * std::log(p) + (1.0 - obs.T) * std::log(1.0 - p);
        }

        r.beta = beta;
        r.logLikelihood = ll;
        for (std::size_t i = 0; i < N; ++i) {
            r.se[i] = std::sqrt(XtWXinv[i][i]); // inverse observed Fisher information
            r.zstat[i] = (r.se[i] > 0.0) ? r.beta[i] / r.se[i] : 0.0;
            r.pvalue[i] = 2.0 * (1.0 - normalCDF(std::fabs(r.zstat[i])));
        }
        return r;
    }

    std::vector<Observation> data_;
};

// ---------- demo / usage ----------

static void printFit(const FitResult& r, long n, std::ostream& out = std::cerr) {
    static const char* names[6] = {"Intercept (b0)", "D slope (b1)",
                                    "C offset (b2)", "D*C interaction (b3)",
                                    "D*D curvature (b4)", "D*D*C interaction (b5)"};
    out << "n = " << n
              << ", converged = " << (r.converged ? "yes" : "no")
              << " (" << r.iterations << " iterations)"
              << ", logLik = " << r.logLikelihood << "\n";
    out.setf(std::ios::fixed);
    out.precision(5);
    for (int i = 0; i < 6; ++i) {
        out << "  " << names[i]
                  << ": beta=" << r.beta[i]
                  << "  se=" << r.se[i]
                  << "  z=" << r.zstat[i]
                  << "  p=" << r.pvalue[i]
                  << "\n";
    }
    out << "  -> Linear slope difference   (C=true - C=false) = b3 = " << r.beta[3]
              << ", p = " << r.pvalue[3]
              << (r.pvalue[3] < 0.05 ? "  => significant at alpha=0.05\n"
                                      : "  => not significant at alpha=0.05\n");
    out << "  -> Curvature difference      (C=true - C=false) = b5 = " << r.beta[5]
              << ", p = " << r.pvalue[5]
              << (r.pvalue[5] < 0.05 ? "  => significant at alpha=0.05\n"
                                      : "  => not significant at alpha=0.05\n");
}

static void printLRT(const LikelihoodRatioTest& t, std::ostream& out = std::cerr) {
    out << "  Likelihood-ratio test (full vs. reduced, df=" << t.df << "):\n";
    out << "    llFull=" << t.llFull << "  llReduced=" << t.llReduced
        << "  LR=" << t.lrStatistic << "  p=" << t.pvalue
        << (t.pvalue < 0.05 ? "  => curves differ significantly\n"
                             : "  => no significant difference between curves\n");
}

/*
int main() {
    IncrementalLogisticRegression reg;

    // ---- Example: synthetic data, added ONE POINT AT A TIME ----
    // True model, C=false: logit(P) = -3.0 + 1.0*D + 0.0*D*D
    // True model, C=true : logit(P) = -3.0 + 3.0*D + 0.4*D*D   (different slope AND curvature!)
    // Replace this generator with your real incremental data source --
    // just call reg.addPoint(D, C, T) whenever a new observation arrives.

    std::mt19937 rng(42); // fixed seed for reproducible demo output
    std::uniform_real_distribution<double> dDist(0.2, 3.0);
    std::uniform_real_distribution<double> uDist(0.0, 1.0);

    const double b0 = -3.0, bFalseLin = 1.0, bTrueLin = 3.0;
    const double bFalseQuad = 0.0, bTrueQuad = 0.4;
    const int perGroup = 60;

    std::vector<Observation> pool;
    for (int i = 0; i < perGroup; ++i) {
        double D = dDist(rng);
        double etaF = b0 + bFalseLin * D + bFalseQuad * D * D;
        double pF = sigmoid(etaF);
        int8_t T_F = (uDist(rng) < pF) ? 1 : 0;
        pool.push_back({static_cast<int8_t>(std::lround(D * 10)), 0, T_F});

        D = dDist(rng);
        double etaT = b0 + bTrueLin * D + bTrueQuad * D * D;
        double pT = sigmoid(etaT);
        int8_t T_T = (uDist(rng) < pT) ? 1 : 0;
        pool.push_back({static_cast<int8_t>(std::lround(D * 10)), 1, T_T});
    }
    // shuffle so points arrive in a realistic mixed order
    std::shuffle(pool.begin(), pool.end(), rng);

    for (const auto& obs : pool) {
        reg.addPoint(obs.D, obs.C, obs.T);

        if (reg.count() % 20 == 0) {
            std::cerr << "---- after " << reg.count() << " points ----\n";
            printFit(reg.fit(), reg.count());
        }
    }

    std::cerr << "\n==== FINAL RESULT ====\n";
    printFit(reg.fit(), reg.count());
    printLRT(reg.likelihoodRatioTest());

    return 0;
}
*/
