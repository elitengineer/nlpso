// Copyright (C) 2022 Robert Hoffmann <robert.hoffmann@smail.emt.h-brs.de>
// I'll release this under a license once I decided which.
#pragma once

#include <random>
#include <limits>
#include <functional>
#include <iostream> // Remove after debugging //

struct extrenum_t
{
	double value;
	double *point;
};

struct nlpso_cfg_t
{
	int xdim, ldim;
	int swarmsize_pp, swarmsize_bif;
	int iterations_pp, iterations_bif;
	int period;
	double mu;
	double Cstop_pp, Cstop_bif;
	double c_inertia, c_personal, c_group;
	double *xmin, *xmax, *lmin, *lmax;
	std::function<double(double*, double*)> f;
	std::function<double(double, double*, int)> objective_pp;
	std::function<double(nlpso_cfg_t, extrenum_t, double*)> objective_bif;
};

extrenum_t PSOpp(nlpso_cfg_t cfg, double *lambda);
extrenum_t PSObif(nlpso_cfg_t cfg);