#ifndef ORDINAL_LOGISTIC_REGRESSION
#define ORDINAL_LOGISTIC_REGRESSION

#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

class OnlineOrdinalAdam {
public:
    struct Sample {
        std::vector<double> x; // predictors
        int y;                 // category 0..K-1
    };

    OnlineOrdinalAdam(int K0, int p0,
                      double lr = 0.01,
                      double b1 = 0.9,
                      double b2 = 0.999,
                      double eps = 1e-8,
		      const std::vector<double>& betaInit = std::vector<double>(),
		      double l2 = 0,
		      double p0_ = 1)
        : K(K0), p(p0), L2(l2), P0(p0_),
          delta(K0 - 1, 1.0 / K0),
          //beta(p0, 0.0),
          beta(betaInit),
          m_delta(K0 - 1, 0.0),
          v_delta(K0 - 1, 0.0),
          m_beta(p0, 0.0),
          v_beta(p0, 0.0),
          learning_rate(lr),
          beta1(b1),
          beta2(b2),
          epsilon(eps),
          t(0),
	  G_L2(0),
          G_delta(K - 1, 0.0),
          G_beta(p, 0.0)	{
		  beta.resize(p0, 0.0);
	  }

    // ---------- thresholds ----------
    std::vector<double> get_alpha() const {
        std::vector<double> a(K - 1);
        double acc = 0.0;
        for (int j = 0; j < K - 1; ++j) {
            acc += std::exp(delta[j]);
            a[j] = acc;
        }
        return a;
    }

    // ---------- sigmoid ----------
    static double sigmoid(double z) {
        if (z >= 0.0)
            return 1.0 / (1.0 + std::exp(-z));
        double e = std::exp(z);
        return e / (1.0 + e);
    }

    // ---------- cumulative P(Y <= j) ----------
    double cumulative(int j,
                      const std::vector<double>& x,
                      const std::vector<double>& alpha) const {
        double xb = dot(beta, x);
        return sigmoid(alpha[j] - xb);
    }

    // ---------- category probability ----------
    double prob(const Sample& s) const {
        auto alpha = get_alpha();
        if (s.y == 0)
            return cumulative(0, s.x, alpha);
        if (s.y == K - 1)
            return 1.0 - cumulative(K - 2, s.x, alpha);
        return cumulative(s.y, s.x, alpha) - cumulative(s.y - 1, s.x, alpha);
    }

    // ============================================================
    //  MINI-BATCH ADAM UPDATE
    // ============================================================
    void update_batch(const std::vector<Sample>& batch) {
        std::vector<double> g_delta(K - 1, 0.0);
        std::vector<double> g_beta(p, 0.0);
        double g_L2 = 0;

        // accumulate gradients over batch
        for (const auto& s : batch) {
            accumulate_gradient(s, g_delta, g_beta, g_L2);
        }

        // average gradient
        for (double& v : g_delta) v /= batch.size();
        for (double& v : g_beta)  v /= batch.size();

        // ADAM update
        adam_update(delta, m_delta, v_delta, g_delta);
        adam_update(beta,  m_beta,  v_beta,  g_beta);
    }

    // ============================================================
    //  SINGLE-SAMPLE UPDATE (mini-batch size = 1)
    // ============================================================
    void update(const Sample& s) {
        std::vector<double> g_delta(K - 1, 0.0);
        std::vector<double> g_beta(p, 0.0);
        double g_L2 = 0;

        accumulate_gradient(s, g_delta, g_beta, g_L2);

        adam_update(delta, m_delta, v_delta, g_delta, 0);
        adam_update(beta,  m_beta,  v_beta,  g_beta, g_L2);
    }

    void addData(const Sample& s) {
        accumulate_gradient(s, G_delta, G_beta, G_L2);
    }

    void update() {
        adam_update(delta, m_delta, v_delta, G_delta, 0);
        adam_update(beta,  m_beta,  v_beta,  G_beta, G_L2);
    }

    std::vector<double> getParams() const
    {
	    std::vector<double> pa = beta;
	    pa.insert(pa.end(), delta.begin(), delta.end());
	    return pa;
    }

    // ---------- accessors ----------
    std::vector<double> get_beta() const { return beta; }

    double getLoss() const { return loss/n; }

    void endIteration() {
        G_delta = std::vector<double>(K - 1, 0.0);
        G_beta = std::vector<double>(p, 0.0);
	G_L2 = 0;
	loss = 0;
	n = 0;
    }

private:
    int K, p;

    double loss = 0;
    double L2 = 0;
    double P0 = 0;

    std::vector<double> delta; // threshold increments
    std::vector<double> beta;  // slopes

    // ADAM state
    std::vector<double> m_delta, v_delta;
    std::vector<double> m_beta,  v_beta;

    double learning_rate;
    double beta1, beta2, epsilon;
    long long t;
    int64_t n = 0;

    double G_L2;
    std::vector<double> G_delta;
    std::vector<double> G_beta;


    // ---------- accumulate gradient for one sample ----------
    void accumulate_gradient(const Sample& s,
                             std::vector<double>& g_delta,
                             std::vector<double>& g_beta,
			     double& g_L2) {
	n++;

        auto alpha = get_alpha();

        std::vector<double> c(K - 1), dcdeta(K - 1);
        double xb = dot(beta, s.x);

        for (int j = 0; j < K - 1; ++j) {
            double e = alpha[j] - xb;
            double sig = sigmoid(e);
            c[j] = sig;
            dcdeta[j] = sig * (1.0 - sig);
        }

        double pk = prob(s);
        if (pk < 1e-12) pk = 1e-12;

	// ADD LOSS ACCUMULATION
	loss += -std::log(pk);
	
        std::vector<double> dP_deta(K - 1, 0.0);
        if (s.y == 0) {
            dP_deta[0] = dcdeta[0];
        } else if (s.y == K - 1) {
            dP_deta[K - 2] = -dcdeta[K - 2];
        } else {
            dP_deta[s.y]     += dcdeta[s.y];
            dP_deta[s.y - 1] -= dcdeta[s.y - 1];
        }

        // thresholds
        for (int j = 0; j < K - 1; ++j) {
            double d_alpha_d_delta = std::exp(delta[j]);
            //g_delta[j] += (dP_deta[j] * d_alpha_d_delta) / pk;
            g_delta[j] -= (dP_deta[j] * d_alpha_d_delta) / pk;
        }

        // slopes
        for (int i = 0; i < p; ++i) {
            double grad = 0.0;
            for (int j = 0; j < K - 1; ++j)
                grad += dP_deta[j] * (-s.x[i]);
            //g_beta[i] += grad / pk;
            g_beta[i] -= grad / pk;
        }

	// regularization
	//double P0 = 1.0; // param start value
	double beta_sum = -p * P0; 
        for (int i = 0; i < p; ++i)
	    beta_sum += beta[i];

	loss += L2 * beta_sum*beta_sum;
	g_L2 += L2 * beta_sum;

	/*
        for (int i = 0; i < p; ++i) {
		loss += L2*beta[i]*beta[i];
		g_L2 += beta[i];
	}
	*/
    }

    // ---------- ADAM update ----------
    void adam_update(std::vector<double>& param,
                     std::vector<double>& m,
                     std::vector<double>& v,
                     const std::vector<double>& g,
		     double g_l2 = 0) {
        t++;

        for (size_t i = 0; i < param.size(); ++i) {
	    double gi = g[i] + g_l2;
            m[i] = beta1 * m[i] + (1 - beta1) * gi;
            v[i] = beta2 * v[i] + (1 - beta2) * gi * gi;

            double m_hat = m[i] / (1 - std::pow(beta1, t));
            double v_hat = v[i] / (1 - std::pow(beta2, t));

            param[i] -= learning_rate * m_hat / (std::sqrt(v_hat) + epsilon);
            //param[i] -= learning_rate * m_hat / (std::sqrt(v_hat) + epsilon) + L2 * g_l2;
        }
    }

    static double dot(const std::vector<double>& a,
                      const std::vector<double>& b) {
        double s = 0.0;
        for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
        return s;
    }
};

#endif
