#include <iostream>
#include <cmath>
#include <limits>
#include <new>
#include <chrono>
#include <random>
#include <functional>
#include "nlpso.hpp"

double f(double *x, double *lambda);
double Fpp(double z0, double *lambda, int period);
double Fbif(double z0, double *lambda, double mu, int period, double Fpp_g, double Cpp);

int main()
{
	//Known Test Parameters
	double x = 0.501492;
	double l[] = {0.50511144, 1.4159996};
	double Cpp = 1e-5;
	//F_pp  = 2.77618e-07
	//F_bif = 9.97771e-05
	double result, Fpp_g;
	Fpp_g = Fpp(x, l, 2);

	// Fbif timing test
	auto start = std::chrono::high_resolution_clock::now();
	result = Fbif(x, l, -1.0, 2, Fpp_g, Cpp);
	auto stop = std::chrono::high_resolution_clock::now();

	std::cout << result << std::endl;
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() << std::endl;
	// Fbif timing test stop

	// Test for PSOpp
	double xxmin = 0.00;
	double xxmax = 0.01;
	nlpso_cfg_t cfg =
	{
		.xdim = 1,
		.particles_pp = 16,
		.xmin = &xxmin,
		.xmax = &xxmax,
		.f = f,
		.weight_pp = Fpp,
		.weight_bif = Fbif
	};
	PSOpp(cfg, l);
	// End of PSOpp test

	return 0;
}

// Discrete-time dynamical systems always must use this form.
// These functions get passed to the NLPSO, that's why they need to have the same input/outputs.
double f(double *x, double *lambda)
{
	return fmod((x[0] + lambda[0] - lambda[1]/(2 * M_PI) * sin(2 * M_PI * x[0])), 1.0);
}

double Fpp(double z0, double *lambda, int period)
{
	double z = z0;
	for (int i = 0; i < period; i++)
	{
		z = f(&z, lambda);
	}
	return fabs(z - z0);
}

double Fbif(double z0, double *lambda, double mu, int period, double Fpp_g, double Cpp)
{
	if (Fpp_g < Cpp)
	{
		double z[period];
		z[0] = z0;
		for (int k = 0; k < (period - 1); k++)
		{
			z[k+1] = f(&z[k], lambda);
		}
		double tmp = 1.0;
		for (int k = 0; k < period; k++)
		{
			tmp = tmp * (1.0 - lambda[1] * cos(2 * M_PI * z[k]));
		}
		return fabs(tmp - mu);
	}
	else
	{
		return std::numeric_limits<double>::max();
	}
}