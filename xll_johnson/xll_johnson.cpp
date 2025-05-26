// xll_johnson.cpp
#include "xll_johnson.h"

using namespace xll;

AddIn xai_calculate_now(
	Macro("xll_calculate_now", "CALCULATE.NOW")
);
int WINAPI xll_calculate_now()
{
#pragma XLLEXPORT
	static OPER eq("=");

	return Excel(xlcFormulaReplace, eq, eq) == true;
}
Auto<Open> xao_calculate_now([]() {
	return Excel(xlcOnSheet, OPER(), OPER("CALCULATE.NOW"), true) == true;
});

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
		Arg(XLL_DOUBLE, "xi", "is the xi parameter.", 0),
		Arg(XLL_DOUBLE, "lambda", "is the lambda parameter.", 1),
		Arg(XLL_DOUBLE, "gamma", "is the gamma parameter.", 0),
		Arg(XLL_DOUBLE, "delta", "is the delta parameter.", 1),
		})
		.Uncalced()
	.Category("XLL")
	.FunctionHelp("Returns a handle to a Johnson distribution.")
);
HANDLEX WINAPI xll_variate_johnson(double xi, double lambda, double gamma, double delta)
{
#pragma XLLEXPORT
	HANDLEX h = INVALID_HANDLEX;

	try {
		handle<fms::Variate> h_(new fms::Johnson(xi, lambda, gamma, delta));
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
		Arg(XLL_DOUBLE, "s", "is the vol."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the cumulative distribution function .")
);
double WINAPI xll_variate_cdf(HANDLEX h, double x, double s)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);

		return h_->cdf(x, s);
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
		Arg(XLL_DOUBLE, "s", "is the vol."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the probability density function of the variate distribution.")
);
double WINAPI xll_variate_pdf(HANDLEX h, double x, double s)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);

		return h_->pdf(x, s);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_variate_cgf(
	Function(XLL_DOUBLE, "xll_variate_cgf", "VARIATE.CGF")
	.Arguments({
		Arg(XLL_HANDLEX, "h", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "s", "is the value at which to evaluate the cgf."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the cumulant generating function of the variate distribution.")
);
double WINAPI xll_variate_cgf(HANDLEX h, double s)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);

		return h_->cgf(s);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_variate_mgf(
	Function(XLL_DOUBLE, "xll_variate_mgf", "VARIATE.MGF")
	.Arguments({
		Arg(XLL_HANDLEX, "h", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "s", "is the value at which to evaluate the mgf."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the moment generating function of the variate distribution.")
);
double WINAPI xll_variate_mgf(HANDLEX h, double s)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);

		return h_->mgf(s);
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

		return j->x(z);
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

AddIn xai_normal_moneyness(
	Function(XLL_DOUBLE, "xll_normal_moneyness", "NORMAL.MONEYNESS")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward price."),
		Arg(XLL_DOUBLE, "s", "is the vol."),
		Arg(XLL_DOUBLE, "k", "is the strike price."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the moneyness of a lognormal distribution.")
);
double WINAPI xll_normal_moneyness(double f, double s, double k)
{
#pragma XLLEXPORT
	try {
		return fms::moneyness(f, s, k);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_normal_put_value(
	Function(XLL_DOUBLE, "xll_normal_put_value", "NORMAL.PUT.VALUE")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward price."),
		Arg(XLL_DOUBLE, "s", "is the vol."),
		Arg(XLL_DOUBLE, "k", "is the strike price."),
		})
	.Category("XLL")
	.FunctionHelp("Returns the put value of a normal distribution.")
);
double WINAPI xll_normal_put_value(double f, double s, double k)
{
#pragma XLLEXPORT
	try {
		return fms::put_value(f, s, k);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_normal_put_delta(
	Function(XLL_DOUBLE, "xll_normal_put_delta", "NORMAL.PUT.DELTA")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward price."),
		Arg(XLL_DOUBLE, "s", "is the vol."),
		Arg(XLL_DOUBLE, "k", "is the strike price."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the put delta of a normal distribution.")
);
double WINAPI xll_normal_put_delta(double f, double s, double k)
{
#pragma XLLEXPORT
	try {
		return fms::put_delta(f, s, k);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_normal_put_implied(
	Function(XLL_DOUBLE, "xll_normal_put_implied", "NORMAL.PUT.IMPLIED")
	.Arguments({
		Arg(XLL_DOUBLE, "f", "is the forward price."),
		Arg(XLL_DOUBLE, "p", "is the put value."),
		Arg(XLL_DOUBLE, "k", "is the strike price."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the put implied of a normal distribution.")
);
double WINAPI xll_normal_put_implied(double f, double p, double k)
{
#pragma XLLEXPORT
	try {
		return fms::put_implied(f, p, k);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	return std::numeric_limits<double>::quiet_NaN();
}

AddIn xai_johnson_put_value(
	Function(XLL_DOUBLE, "xll_johnson_put_value", "JOHNSON.PUT.VALUE")
	.Arguments({
		Arg(XLL_HANDLEX, "h", "is the handle to a Johnson distribution."),
		Arg(XLL_DOUBLE, "k", "is the strike price."),
		})
		.Category("XLL")
	.FunctionHelp("Returns the put value of a johnson distribution.")
);
double WINAPI xll_johnson_put_value(HANDLEX h, double k)
{
#pragma XLLEXPORT
	try {
		handle<fms::Variate> h_(h);
		ensure(h_);
		const auto j = h_.as<fms::Johnson>();
		ensure(j);

		return fms::put_value(*j, k);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	return std::numeric_limits<double>::quiet_NaN();
}
