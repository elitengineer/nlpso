// Copyright (C) 2022 Robert Hoffmann <robert.hoffmann@smail.emt.h-brs.de>
// I'll release this under a license once I decided which.
#pragma once

#include <random>
#include <limits>
#include <functional>
#include <array>
#include <cstring>
#include <iostream> // Remove after debugging //
#include <chrono>

class particle
{
private:
	int dim;
public:
	std::vector<double> position;
	std::vector<double> pbest;
	std::vector<double> velocity;
	double Fposition = 0;
	double Fpbest = std::numeric_limits<double>::max();

	// Having the constructor do the new operation would be better, but I don't know
	// how to suppress the creation of a new object when initializing the array
	// which leads to double freeing.
	particle(){}

	void init(const int &dimensions, std::default_random_engine &gen, std::uniform_real_distribution<> searchspace[])
	{
		dim = dimensions;
		position.reserve(dim);
		pbest.reserve(dim);
		velocity.reserve(dim);
		for (int d = 0; d < dim; d++)
		{
			position.push_back(searchspace[d](gen));
		}
		pbest = position;
		std::fill(velocity.begin(), velocity.end(), 0);
	}

	bool is_within_bounds(double *min, double *max)
	{
		for (int d = 0; d < dim; d++)
		{
			if ((position[d] < min[d]) || (position[d] > max[d]))
			{
				return false;
			}
		}
		return true;
	}

	void update_personal()
	{
		if (Fposition < Fpbest)
		{
			Fpbest = Fposition;
			for (int d = 0; d < dim; d++) // Elementwise is faster
			{
				pbest[d] = position[d];
			}
		}
	}

	void update_position(double w, double cp, double cg, double *gbest)
	{
		for (int d = 0; d < dim; d++)
		{
			velocity[d] =	w * velocity[d]
							+ cp * (pbest[d] - position[d])
							+ cg * (gbest[d] - position[d]);

			position[d] = position[d] + velocity[d];
		}
	}
};

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
	std::function<double(double*, const std::vector<double>&)> f;
	std::function<double(const double&, const std::vector<double>&, const int&)> objective_pp;
	std::function<double(nlpso_cfg_t, const extrenum_t&, const std::vector<double>&)> objective_bif;
};

extrenum_t PSOpp(nlpso_cfg_t cfg, const std::vector<double> &lambda);
extrenum_t PSObif(const nlpso_cfg_t &cfg);