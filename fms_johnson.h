// fms_johnson.h - Johnson distribution
#include <cmath>
#include <numbers>

namespace fms {

	struct Variate {
		virtual ~Variate() = default;
		// Cumulative share distribution E[exp(X - kappa(s)) 1(X <= x)].
		double cdf(double x, double s = 0) const
		{
			return _cdf(x, s);
		}
		// Share probability density d/dx cdf(x, s).
		double pdf(double x, double s = 0) const
		{
			return _pdf(x, s);
		}
		// Cumulant generating function log E[exp(s X)].
		double cgf(double s)
		{
			return _cgf(s);
		}
		// Moment generating function E[exp(s X)].
		double mgf(double s)
		{
			return _mgf(s);
		}
	private:
		virtual double _pdf(double x, double s) const = 0;
		virtual double _cdf(double x, double s) const = 0;
		virtual double _cgf(double s) const = 0;
		virtual double _mgf(double s) const = 0;
	};

	struct Normal : public Variate {
		double mu, sigma;

		Normal(double mu = 0, double sigma = 1)
			: mu(mu), sigma(sigma)
		{
			if (sigma <= 0) {
				throw std::invalid_argument("sigma must be positive");
			}
		}
		double _pdf(double x, double s) const override
		{
			double z = (x - s - mu) / sigma;

			return std::exp(-z * z / 2) * std::numbers::inv_sqrtpi / (sigma * std::numbers::sqrt2);
		}
		double _cdf(double x, double s) const override
		{
			double z = (x - s - mu) / sigma;

			return 0.5 * (1 + std::erf(z / std::numbers::sqrt2));
		}
		double _cgf(double s) const override
		{
			return mu * s + 0.5 * sigma * sigma * s * s;
		}
		double _mgf(double s) const override
		{
			return std::exp(_cgf(s));
		}
	};

	// X is Johnson S_U if
	// N = gamma + delta * asinh((X - xi) / lambda) is normal
	class Johnson : public Variate {
#ifdef _DEBUG
	public:
#endif
		double gamma, delta, lambda, xi;
		Normal N;

		// To normal.
		double n(double x) const
		{
			return gamma + delta * std::asinh((x - xi) / lambda);
		}
		// dN/dX
		double dn_dx(double x) const
		{
			double y = (x - xi) / lambda;

			return delta / (lambda * std::sqrt(1 + y * y));
		}
		// From normal
		double x(double n) const
		{
			return xi + lambda * std::sinh((n - gamma) / delta);
		}

		// E[(sinh^k((Z - gamma)/delta)] = 2^{-k} sum_{j=0}^{k} (-1)^j C_kj exp((k - 2j)^2/delta^2) - (k - 2j)gamma/delta)
		double Esinh_k(int k) const
		{
			double E = 0;
			double _C_kj = 1; // (-1)^j k choose j
			int k_ = k;

			for (int j = 0; j <= k; ++j) {
				double k_2j = k - 2 * j;
				E += _C_kj * std::exp(k_2j * k_2j / (2 * delta * delta) - k_2j * gamma / delta);
				_C_kj *= -k_--;
				_C_kj /= j + 1;
			}

			return E * std::pow(2.0, -k);
		}
		// E[X^n] = E[(xi + lambda * sinh((Z - gamma)/delta))^n]
		//        = sum_{k=0}^{n} C_nk xi^(n - k) lambda^k E[sinh^k((Z - gamma)/delta)]
		double _moment(unsigned n) const
		{
			double mn = 0;

			double Cnk = 1;
			double n_ = n;
			for (unsigned k = 0; k <= n; ++k) {
				mn += Cnk * pow(xi, n - k) * pow(lambda, k) * Esinh_k(k);
				Cnk *= n_--;
				Cnk /= k + 1;
			}

			return mn;
		}
	public:
		Johnson(double gamma, double delta, double lambda, double xi, double mu = 0, double sigma = 1)
			: gamma(gamma), delta(delta), lambda(lambda), xi(xi), N(mu, sigma)
		{
			if (delta <= 0) {
				throw std::invalid_argument("delta must be positive");
			}
			if (lambda <= 0) {
				throw std::invalid_argument("lambda must be positive");
			}
		}
		// E[X]
		double mean() const
		{
			return xi + lambda * std::exp(N.sigma * N.sigma / (2 * delta * delta)) * std::sinh((N.mu - gamma) / delta);
		}
		// set E[X] to f
		void mean(double f)
		{
			xi += f - mean();
		}
		// Probability density function.
		double _pdf(double x, double s = 0) const override
		{
			return N.pdf(n(x)) * dn_dx(x);
		}
		// Cumulative distribution function.
		double _cdf(double x, double s = 0) const override
		{
			return N.cdf(n(x));
		}
		double _cgf(double s) const override
		{
			return std::log(_mgf(s));
		}
		double _mgf(double s) const override
		{
			double mgf = 1; // first term

			double sn = 1; // s^n/n!
			int dup = 1; // protect against 0 moments
			for (unsigned n = 1; n < 100; ++n) {
				sn *= s / n;
				double incr = moment(n) * sn;
				mgf += incr;
				if (abs(incr) < 1e-8) {
					if (dup == 0) {
						break;
					}
					else {
						--dup;
					}
				}
			}

			return mgf;
		}

		double moment(unsigned k) const
		{
			return _moment(k);
		}

		// Share cumulative distribution E[X 1(X <= x)].
		double cdf_(double x) const
		{
			double z = (N.mu - gamma) / delta;
			double s = N.sigma / delta;
			double dN = std::exp(z) * N.cdf(n(x) - s*s * delta) - std::exp(-z) * N.cdf(n(x) + s*s * delta);

			return xi * N.cdf(n(x)) + 0.5 * lambda * std::exp(s * s / 2) * dN ;
		}

	};

	double put_value(const Johnson& j, double k)
	{
		return k * j.cdf(k) - j.cdf_(k);
	}
	// B-S/M
	double moneyness(double f, double s, double k)
	{
		if (s <= 0) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return (std::log(k / f) + s * s / 2) / s;
	}

	double put_value(double f, double s, double k)
	{
		static Normal N;
		double z = moneyness(f, s, k);

		return k * N.cdf(z) - f * N.cdf(z - s);
	}
	double put_delta(double f, double s, double k)
	{
		static Normal N;
		double z = moneyness(f, s, k);

		return -N.cdf(z, s);
	}
	double put_vega(double f, double s, double k)
	{
		static Normal N;
		double z = moneyness(f, s, k);

		return f * N.pdf(z, s);
	}

	double put_implied(double f, double p, double k, double s = 0.1, double eps = 1e-8, size_t iter = 100)
	{
		while (iter--) {
			double s_ = s - (put_value(f, s, k) - p) / put_vega(f, s, k);
			if (s_ < 0) {
				s_ = s/2;
			}
			if (std::fabs(s_ - s) < eps) {
				break;
			}
			s = s_;
		}

		return s;
	}

}