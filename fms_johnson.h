// fms_johnson.h - Johnson distribution
#include <cmath>
#include <numbers>

namespace fms {

	struct Variate {
		double pdf(double x) const
		{
			return _pdf(x);
		}
		double cdf(double x) const
		{
			return _cdf(x);
		}
	private:
		virtual double _pdf(double) const = 0;
		virtual double _cdf(double) const = 0;	
	};

	class Normal : public Variate {
		double mu, sigma;
	public:
		Normal(double mu = 0, double sigma = 1)
			: mu(mu), sigma(sigma) 
		{ 
			if (sigma <= 0) {
				throw std::invalid_argument("sigma must be positive");
			}
		}
		double _pdf(double x) const override
		{
			double z = (x - mu) / sigma;

			return std::exp(-z*z/2) * std::numbers::inv_sqrtpi / std::numbers::sqrt2;;
		}
		double _cdf(double x) const override
		{
			double z = (x - mu) / sigma;

			return 0.5 * (1 + std::erf(z / std::numbers::sqrt2));
		}
	};

	// X is Johnson S_U if
	// Z = gamma + delta * asinh((X - xi) / lambda) is normal
	class Johnson : public Variate {
#ifdef _DEBUG
	public:
#endif

		double gamma, delta, lambda, xi;

		// To std normal.
		double Z(double x) const
		{
			return gamma + delta * std::asinh((x - xi) / lambda);
		}
		// dZ/dX
		double dZdX(double x) const
		{
			double y = (x - xi) / lambda;

			return delta / (lambda * std::sqrt(1 + y * y));
		}
		// From std normal
		double X(double z) const
		{
			return xi + lambda * std::sinh((z - gamma) / delta);
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
		Johnson(double gamma, double delta, double lambda, double xi)
			: gamma(gamma), delta(delta), lambda(lambda), xi(xi) 
		{
			if (delta <= 0) {
				throw std::invalid_argument("delta must be positive");
			}
			if (lambda <= 0) {
				throw std::invalid_argument("lambda must be positive");
			}
		}
		// Probability density function.
		double _pdf(double x) const override
		{
			return Normal(0,1).pdf(Z(x)) * dZdX(x);
		}
		// Cumulative distribution function.
		double _cdf(double x) const override
		{
			return Normal(0, 1).cdf(Z(x));
		}
		double moment(unsigned k) const
		{
			return _moment(k);
		}
	};
}