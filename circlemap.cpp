// Copyright (C) 2022 Robert Hoffmann <robert.hoffmann@smail.emt.h-brs.de>
// I'll release this under a license once I decided which.
#include <iostream>
#include <cmath>
#include <limits>
#include <new>
#include <chrono>
#include <random>
#include <functional>
#include <fstream>
#include "nlpso.hpp"

double f(double *x, double *lambda);
double Fpp(double z0, double *lambda, int period);
double Fbif(nlpso_cfg_t cfg, extrenum_t pp, double *lambda);

int main()
{
	//Known Test Parameters
	//double x = 0.501492;
	double l[] = {0.50511144, 1.4159996};
	//F_pp  = 2.77618e-07
	//F_bif = 9.97771e-05

	// Test for PSObif and PSOpp
	double xxmin = 0.00;
	double xxmax = 1.00;
	double llmin[2] = {0.35, 0.50};
	double llmax[2] = {0.65, 2.00};
	nlpso_cfg_t cfg =
	{
		.xdim = 1,
		.ldim = 2,
		.swarmsize_pp = 30,
		.swarmsize_bif = 30,
		.iterations_pp = 300,
		.iterations_bif = 300,
		.period = 2,
		.mu = -1.0,
		.Cstop_pp = 1e-5,
		.Cstop_bif = 1e-3,
		.c_inertia = 0.729,
		.c_personal = 1.494,
		.c_group = 1.494,
		.xmin = &xxmin,
		.xmax = &xxmax,
		.lmin = llmin,
		.lmax = llmax,
		.f = f,
		.objective_pp = Fpp,
		.objective_bif = Fbif
	};
	// Test for PSOpp
	double test = PSOpp(cfg, l).point[0];
	std::cout << test << std::endl;
	std::cout << Fpp(test, l, 2) << std::endl;
	// End of PSOpp test
	
	// Test for PSObif
	int points = 100;
	extrenum_t xp[points];
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < points; i++)
	{
		xp[i] = PSObif(cfg);
	}
	auto stop = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << std::endl;
	
	// End of PSObif test

	// Writing to CSV file for plotting in MATLAB
	
	std::ofstream plotpoints;
	plotpoints.open("plotpoints.csv");
	for (int d = 0; d < cfg.ldim; d++)
	{
		for (int i = 0; i < points; i++)
		{
			plotpoints << xp[i].point[d] << ",";
		}
		plotpoints << "\n";
	}
	plotpoints.close();
	
	// End of CSV business
	
	// Yes, I am not deleting[] the xp.point values, since this program ends here anyway.
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

double Fbif(nlpso_cfg_t cfg, extrenum_t pp, double *lambda)
{
	if (pp.value < cfg.Cstop_pp)
	{
		double z[cfg.period];
		z[0] = pp.point[0];
		for (int k = 0; k < (cfg.period - 1); k++)
		{
			z[k+1] = f(&z[k], lambda);
		}
		double tmp = 1.0;
		for (int k = 0; k < cfg.period; k++)
		{
			tmp = tmp * (1.0 - lambda[1] * cos(2 * M_PI * z[k]));
		}
		return fabs(tmp - cfg.mu);
	}
	else
	{
		return std::numeric_limits<double>::max();
	}
}