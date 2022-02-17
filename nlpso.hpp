// Copyright (C) 2022 Robert Hoffmann <robert.hoffmann@smail.emt.h-brs.de>
// I'll release this under a license once I decided which.
#pragma once

#include <random>
#include <limits>
#include <functional>
#include <iostream> // Remove after debugging //

struct nlpso_cfg_t
{
	int xdim, ldim;
	int swarmsize_pp, swarmsize_bif;
	int iterations_pp, iterations_bif;
	int period;
	double Cstop_pp, Cstop_bif;
	double c_inertia, c_personal, c_group;
	double *xmin, *xmax, *lmin, *lmax;
	std::function<double(double*, double*)> f;
	std::function<double(double, double*, int)> objective_pp;
	std::function<double(double, double*, double, int, double, double)> objective_bif;
};

struct extrenum_t
{
	double value;
	double *point;
};

extrenum_t PSOpp(nlpso_cfg_t cfg, double *lambda);