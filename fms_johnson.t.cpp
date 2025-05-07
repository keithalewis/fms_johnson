// fms_johnson.t.cpp - Test for Johnson distribution
#include <cassert>
#include <functional>
#include "fms_johnson.h"

double symmetric_difference(const std::function<double(double)>& f, double x, double h)
{
	return (f(x + h) - f(x - h)) / (2 * h);
}

int dist_test()
{
	fms::Johnson j(0, 1, 1, 0);
	{
		double data[][2] = {
			// x, h
			{0.0, 1e-8},
			{1.0, 1e-8},
			{-1.0, 1e-8},
		};
		for (const auto [x, h] : data) {
			double pdf = j.pdf(x);
			double pdf_ = symmetric_difference([&j](double x) { return j.cdf(x); }, x, h);
			double err = pdf - pdf_;
			assert(std::fabs(err) < h);
		}

	}

	return 0;
}

int main()
{
	int test_dist = dist_test();

	return 0;
}