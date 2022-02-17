#include <random>
#include <limits>
#include <functional>
#include <iostream>

struct nlpso_cfg_t
{
	int xdim, ldim;
	int swarmsize_pp, swarmsize_bif;
	int particles_pp, particles_bif;
	int iterations_pp, iterations_bif;
	int period;
	double Cstop_pp, Cstop_bif;
	double c_inertia, c_personal, c_group;
	double *xmin, *xmax, *lmin, *lmax;
	std::function<double(double*, double*)> f;
	std::function<double(double, double*, int)> weight_pp;
	std::function<double(double, double*, double, int, double, double)> weight_bif;
};

double PSOpp(nlpso_cfg_t cfg, double *lambda);