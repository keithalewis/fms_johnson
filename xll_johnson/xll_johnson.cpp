// xll_johnson.cpp
#include "xll_johnson.h"

using namespace xll;

AddIn xai_johnson(
	Function(XLL_HANDLEX, "xll_johnson", "JOHNSON")
	.Arguments({
		Arg(XLL_DOUBLE, "gamma", "is the gamma parameter."),
		Arg(XLL_DOUBLE, "delta", "is the delta parameter."),
		Arg(XLL_DOUBLE, "lambda", "is the lambda parameter."),
		Arg(XLL_DOUBLE, "xi", "is the xi parameter."),
	})
	.Uncalced()
	.Category("XLL")
	.FunctionHelp("Returns the Johnson distribution.")
);
HANDLEX WINAPI xll_johnson(double gamma, double delta, double lambda, double xi)
{
#pragma XLLEXPORT
	HANDLEX h = INVALID_HANDLEX;

	try {
		handle<fms::Johnson> j(new fms::Johnson(gamma, delta, lambda, xi));
		h = j.get();
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}

	return h;
}