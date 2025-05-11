// fms_johnson.t.cpp - Test for Johnson distribution
#include <cassert>
#include <functional>
#include <random>
#include "fms_johnson.h"

static std::default_random_engine e;
static std::uniform_real_distribution<double> u;
static std::normal_distribution<double> n(0, 1);

using fun = std::function<double(double)>;

// (f(x + h) - f(x - h)) / (2 * h) = f'(x) + f'''(x) h^2/3! + O(h^3)
double symmetric_difference(const fun& f, double x, double h)
{
	return (f(x + h) - f(x - h)) / (2 * h);
}

// Third derivative at x.
double bigO(const fun& df, const fun& f, double x, double h)
{
	return (df(x) - symmetric_difference(f, x, h)) / (h * h / 6);
}

// Mean and standard deviation of calling f() N times.
std::pair<double, double> monte(const std::function<double(void)>& f, int N)
{
	double m1 = 0;
	double m2 = 0;

	for (int n = 1; n <= N; ++n) {
		double x = f();
		m1 += (x - m1) / n;
		m2 += (x * x - m2) / n;
	}

	return { m1, sqrt(m2 - m1 * m1) };
}
auto monte_uniform(int n)
{
	const auto f = []() { return u(e); };
	return monte(f, n);
}
int monte_uniform_test(int n)
{
	{
		auto [m, s] = monte_uniform(n);
		assert(std::fabs(m - 0.5) <= 1. / std::sqrt(n));
		assert(std::fabs(s - 1. / sqrt(12)) <= 1. / std::sqrt(n));
	}
	return 0;
}

int normal_test()
{
	{
		fms::Normal j(0, 1);
		double data[][3] = {
			// x, h, O
			{0.0, 1e-4, 1},
			{1.0, 1e-4, 1},
			{-1.0, 1e-4, 1},
		};
		const auto jpdf = [&j](double x) { return j.pdf(x); };
		const auto jcdf = [&j](double x) { return j.cdf(x); };
		for (const auto [x, h, O] : data) {
			double O_ = bigO(jpdf, jcdf, x, h);
			assert(std::fabs(O_) < O);
		}
	}

	return 0;
}

int dist_test()
{
	fms::Johnson j(0, 1, 1, 0);
	{
		double data[][3] = {
			// x, h, O
			{0.0, 1e-4, 1},
			{1.0, 1e-4, 1},
			{-1.0, 1e-4, 1},
		};
		const auto jpdf = [&j](double x) { return j.pdf(x); };
		const auto jcdf = [&j](double x) { return j.cdf(x); };
		for (const auto [x, h, O] : data) {
			double O_ = bigO(jpdf, jcdf, x, h);
			assert(std::fabs(O_) < O);
		}

	}

	return 0;
}

int Esinh_k_test()
{
	double data[][4] = {
		// xi, lambda, gamma, delta
		{0, 1, 1, 0},
		//{1, 1, 1, 0},
		//{0, 1, 1, 1},
		//{1, 1, 1, 1},
	};
	for (auto [xi, lambda, gamma, delta] : data) {
		fms::Johnson j(xi, lambda, gamma, delta);
		int N = 1000000;
		double eps = 2. / sqrt(N);
		static std::normal_distribution<double> n(-j.gamma / j.delta, 1 / j.delta);
		{
			double j1 = j.Esinh_k(1);
			const auto f = []() { return sinh(n(e)); };
			auto [j1_, s] = monte(f, N);
			double err = (j1 - j1_) / s;
			assert(fabs(err) < eps);
		}
		{
			double j2 = j.Esinh_k(2);
			const auto f = []() { return pow(sinh(n(e)), 2); };
			auto [j2_, s] = monte(f, N);
			double err = (j2 - j2_) / s;
			assert(fabs(err) < eps);
		}
		{
			double j3 = j.Esinh_k(3);
			const auto f = []() { return pow(sinh(n(e)), 3); };
			auto [j3_, s] = monte(f, N);
			double err = (j3 - j3_) / s;
			assert(fabs(err) < eps);
		}
		{
			double j4 = j.Esinh_k(4);
			const auto f = []() { return pow(sinh(n(e)), 4); };
			auto [j4_, s] = monte(f, N);
			double err = (j4 - j4_) / s;
			assert(fabs(err) < eps);
		}
	}
	return 0;
}

int main()
{
	assert(!monte_uniform_test(10000));
	assert(!normal_test());
	assert(!dist_test());
	assert(!Esinh_k_test());

	return 0;
}