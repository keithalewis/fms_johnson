// xll_johnson.cpp
#include "xll_johnson.h"

using namespace xll;

AddIn xai_variate_normal(
	Function(XLL_HANDLEX, "xll_variate_normal", "\\VARIATE.NORMAL")
	.Arguments({
		Arg(XLL_DOUBLE, "mu", "is mean.", 0),
		Arg(XLL_DOUBLE, "sigma", "is the standard deviation.", 1),
		})
		.Uncalced()
	.Category("XLL")
	.FunctionHelp("Returns a handle to a normal distribution.")
);
HANDLEX WINAPI xll_variate_normal(double mu, double sigma)
{
#pragma XLLEXPORT
	HANDLEX h = INVALID_HANDLEX;

	try {
		if (sigma <= 0) {
			sigma = 1;
		}
		handle<fms::Variate> h_(new fms::Normal(mu, sigma));
		h = h_.get();
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return h;
}

AddIn xai_variate_johnson(
	Function(XLL_HANDLEX, "xll_variate_johnson", "\\VARIATE.JOHNSON")
	.Arguments({
		Arg(XLL_DOUBLE, "gamma", "is the gamma parameter."),
		Arg(XLL_DOUBLE, "delta", "is the delta parameter."),
		Arg(XLL_DOUBLE, "lambda", "is the lambda parameter."),
		Arg(XLL_DOUBLE, "xi", "is the xi parameter."),
		})
		.Uncalced()
	.Category("XLL")
	.FunctionHelp("Returns a handle to a Johnson distribution.")
);
HANDLEX WINAPI xll_variate_johnson(double gamma, double delta, double lambda, double xi)
{
#pragma XLLEXPORT
	HANDLEX h = INVALID_HANDLEX;

	try {
		handle<fms::Variate> h_(new fms::Johnson(gamma, delta, lambda, xi));
		h = h_.get();
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return h;
}

AddIn xai_variate_cdf(
	Function(XLL_DOUBLE, "xll_variate_cdf", "VARIATE.CDF")
	.Arguments({
		Arg(XLL_HANDLEX, "h", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "x", "is the value at which to evaluate the CDF."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the cumulative distribution function .")
);
double WINAPI xll_variate_cdf(HANDLEX h, double x)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);

		return h_->cdf(x);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_variate_pdf(
	Function(XLL_DOUBLE, "xll_variate_pdf", "VARIATE.PDF")
	.Arguments({
		Arg(XLL_HANDLEX, "h", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "x", "is the value at which to evaluate the PDF."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the probability density function of the variatel distribution.")
);
double WINAPI xll_variate_pdf(HANDLEX h, double x)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);

		return h_->pdf(x);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_johnson_X(
	Function(XLL_DOUBLE, "xll_johnson_X", "JOHNSON.X")
	.Arguments({
		Arg(XLL_HANDLEX, "h", "is a handle to a Johnson distribution."),
		Arg(XLL_DOUBLE, "z", "is value of the underlying normal"),
		})
		.Category("XLL")
	.FunctionHelp("Returns the Johnson distribution at z.")
);
double WINAPI xll_johnson_X(HANDLEX h, double z)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);
		const auto j = h_.as<fms::Johnson>();
		ensure(j);

		return j->X(z);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_johnson_Esinh_k(
	Function(XLL_DOUBLE, "xll_johnson_Esinh_k", "JOHNSON.Esinh_k")
	.Arguments({
		Arg(XLL_HANDLEX, "h", "is a handle to a Johnson distribution."),
		Arg(XLL_UINT, "n", "is expected value of sinh^k(X)."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the expected value of sinh^k(X) for the Johnson distribution.")
);
double WINAPI xll_johnson_Esinh_k(HANDLEX h, unsigned k)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);
		const auto j = h_.as<fms::Johnson>();
		ensure(j);

		return j->Esinh_k(k);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_johnson_moment(
	Function(XLL_DOUBLE, "xll_johnson_moment", "JOHNSON.MOMENT")
	.Arguments({
		Arg(XLL_HANDLEX, "h", "is a handle to a Johnson distribution."),
		Arg(XLL_UINT, "n", "is n-th moment."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the n-th moment of the Johnson distribution.")
);

double WINAPI xll_johnson_moment(HANDLEX h, unsigned n)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);
		const auto j = h_.as<fms::Johnson>();
		ensure(j);

		return j->moment(n);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return std::numeric_limits<double>::quiet_NaN();
}