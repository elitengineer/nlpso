// Copyright (C) 2022 Robert Hoffmann <robert.hoffmann@smail.emt.h-brs.de>
// I'll release this under a license once I decided which.

// Both PSO functions are almost the same,
// only the PSObif has additionally the PSOpp inside it.
// Maybe those 2 functions can be combined into one,
// but this comes later.
#include "nlpso.hpp"

extrenum_t PSOpp(nlpso_cfg_t cfg, double *lambda)
{
	// Seeding rand function. I just copied this from cppreference lol
	std::random_device rd;
	std::mt19937 gen(rd());
	// Start of Initialization
	int index_argmin;
	int iterations = 0;
	double Fgbest = std::numeric_limits<double>::max();
	double *gbest = new double[cfg.xdim]; // Maybe a shared pointer in the future?
	// Assign random position one by one within bounds. May need optimization.
	std::uniform_real_distribution<> searchspace[cfg.xdim];
	for (int d = 0; d < cfg.xdim; d++)
	{
		searchspace[d] = std::uniform_real_distribution<>(cfg.xmin[d], cfg.xmax[d]);
	}
	particle swarm[cfg.swarmsize_pp];
	for (int p = 0; p < cfg.swarmsize_pp; p++)
	{
		swarm[p].init(cfg.xdim, gen, searchspace);
	}
	// End of Initialization

	std::uniform_real_distribution<> variation(0.0, 1.0); // Used later in next position calculation
	while ((iterations < cfg.iterations_pp) && (Fgbest > cfg.Cstop_pp))
	{
		iterations++;
		// Check for out of bound positions and don't apply Fpp to them
		for (int p = 0; p < cfg.swarmsize_pp; p++)
		{
			if(swarm[p].is_within_bounds(cfg.xmin, cfg.xmax))
			{
				swarm[p].Fposition = cfg.objective_pp(swarm[p].position[0], lambda, cfg.period);
			}
		}
		// Check if value of position is lower than pbest best and overwrite them if yes
		for (int p = 0; p < cfg.swarmsize_pp; p++)
		{
			swarm[p].update_personal(); // Can be integrated into upper loop
		}
		// Find index of lowest pbest
		index_argmin = 0;
		double min = swarm[0].Fpbest;
		for (int p = 1; p < cfg.swarmsize_pp; p++)
		{
			if (swarm[p].Fpbest < min)
			{
				min = swarm[p].Fpbest;
				index_argmin = p;
			}
		}
		// Assign gbest from best pbest
		Fgbest = swarm[index_argmin].Fpbest;
		for (int d = 0; d < cfg.xdim; d++)
		{
			gbest[d] = swarm[index_argmin].pbest[d];
		}
		// Calculate velocities and new positions combined with random variation
		for (int p = 0; p < cfg.swarmsize_pp; p++)
		{
			swarm[p].update_position(cfg.c_inertia, cfg.c_personal * variation(gen), cfg.c_group * variation(gen), gbest);
		}
	}

	extrenum_t best
	{
		.value = Fgbest,
		.point = gbest
	};

	return best; // Don't forget delete[] outside this call! (:
}

extrenum_t PSObif(nlpso_cfg_t cfg)
{
	// Seeding rand function. I just copied this from cppreference lol
	std::random_device rd;
	std::mt19937 gen(rd());
	// Start of Initialization
	int index_argmin;
	int iterations = 0;
	double Fgbest = std::numeric_limits<double>::max();
	double *gbest = new double[cfg.ldim]; // Maybe a shared pointer in the future?
	// Assign random position one by one within bounds. May need optimization.
	std::uniform_real_distribution<> searchspace[cfg.ldim];
	for (int d = 0; d < cfg.ldim; d++)
	{
		searchspace[d] = std::uniform_real_distribution<>(cfg.lmin[d], cfg.lmax[d]);
	}
	particle swarm[cfg.swarmsize_bif];
	for (int p = 0; p < cfg.swarmsize_bif; p++)
	{
		swarm[p].init(cfg.ldim, gen, searchspace);
	}
	// End of Initialization

	std::uniform_real_distribution<> variation(0.0, 1.0); // Used later in next position calculation
	while ((iterations < cfg.iterations_bif) && (Fgbest > cfg.Cstop_bif))
	{
		iterations++;
		// Check for out of bound positions and don't apply Fpp to them
		for (int p = 0; p < cfg.swarmsize_bif; p++)
		{
			if(swarm[p].is_within_bounds(cfg.lmin, cfg.lmax))
			{
				extrenum_t xp = PSOpp(cfg, swarm[p].position);
				swarm[p].Fposition = cfg.objective_bif(cfg, xp, swarm[p].position);
				delete[] xp.point; // This memory gets called in PSOpp!
			}
		}
		// Check if value of position is lower than pbest best and overwrite them if yes
		for (int p = 0; p < cfg.swarmsize_bif; p++)
		{
			swarm[p].update_personal(); // Can be integrated into upper loop
		}
		// Find index of lowest pbest
		index_argmin = 0;
		double min = swarm[0].Fpbest;
		for (int p = 1; p < cfg.swarmsize_bif; p++)
		{
			if (swarm[p].Fpbest < min)
			{
				min = swarm[p].Fpbest;
				index_argmin = p;
			}
		}
		// Assign gbest from best pbest
		Fgbest = swarm[index_argmin].Fpbest;
		for (int d = 0; d < cfg.ldim; d++)
		{
			gbest[d] = swarm[index_argmin].pbest[d];
		}
		// Calculate velocities and new positions combined with random variation
		for (int p = 0; p < cfg.swarmsize_bif; p++)
		{
			swarm[p].update_position(cfg.c_inertia, cfg.c_personal * variation(gen), cfg.c_group * variation(gen), gbest);
		}
	}

	extrenum_t best
	{
		.value = Fgbest,
		.point = gbest
	};

	return best; // Don't forget delete[] outside this call! (:
}