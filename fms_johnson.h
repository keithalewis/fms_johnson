// fms_johnson.h - Johnson distribution
#include <cmath>
#include <numbers>

namespace fms {

	// X is Johnson S_U if
	// Z = gamma + delta * asinh((X - xi) / lambda) is normal
	class Johnson {
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
			return delta / (std::sqrt (1 + y * y) * lambda);
		}
		// From std normal
		double X(double z) const
		{
			return xi + lambda * std::sinh((z - gamma) / delta);
		}
		// standard normal CDF
		static double Phi(double z)
		{
			return 0.5 * (1 + std::erf(z / std::numbers::sqrt2));
		}
		// standard normal PDF
		static double phi(double z)
		{
			return std::exp(-0.5 * z * z) * std::numbers::inv_sqrtpi/std::numbers::sqrt2;
		}
	public:
		Johnson(double gamma, double delta, double lambda, double xi)
			: gamma(gamma), delta(delta), lambda(lambda), xi(xi) {
		}
		// Probability density function.
		double pdf(double x) const {
			return phi(Z(x)) * dZdX(x);	
		}
		// Cumulative distribution function.
		double cdf(double x) const {
			return Phi(Z(x));
		}
	};
}