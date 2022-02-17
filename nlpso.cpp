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
	double Fposition[cfg.swarmsize_pp] = {0};
	double Fpbest[cfg.swarmsize_pp];
	for (int p = 0; p < cfg.swarmsize_pp; p++)
	{
		Fpbest[p] = std::numeric_limits<double>::max();
	}
	double Fgbest = std::numeric_limits<double>::max();
	double *gbest = new double[cfg.xdim]; // Maybe a shared pointer in the future?
	double position[cfg.xdim][cfg.swarmsize_pp];
	double pbest[cfg.xdim][cfg.swarmsize_pp];
	double velocity[cfg.xdim][cfg.swarmsize_pp] = {0};
	// Assign random position one by one within bounds. May need optimization.
	for (int d = 0; d < cfg.xdim; d++)
	{
		std::uniform_real_distribution<> searchspace(cfg.xmin[d], cfg.xmax[d]);
		for (int p = 0; p < cfg.swarmsize_pp; p++)
		{
			position[d][p] = searchspace(gen);
		}
		std::copy(&position[d][0], &position[d][cfg.swarmsize_pp], pbest[0]); // Save memcpy z to p
	}
	// End of Initialization

	std::uniform_real_distribution<> variation(0.0, 1.0); // Used later in next position calculation
	while ((iterations < cfg.iterations_pp) && (Fgbest > cfg.Cstop_pp))
	{
		iterations++;
		bool isOutOfBound = false;
		// Check for out of bound positions and don't apply Fpp to them
		for (int p = 0; p < cfg.swarmsize_pp; p++)
		{
			for (int d = 0; d < cfg.xdim; d++)
			{
				if ((position[d][p] < cfg.xmin[d]) || (position[d][p] > cfg.xmax[d]))
				{
					isOutOfBound = true;
					break;
				}
			}
			if (!isOutOfBound)
			{
				isOutOfBound = false;
				Fposition[p] = cfg.objective_pp(position[0][p], lambda, cfg.period);
			}
		}
		// Check if value of position is lower than pbest best and overwrite them if yes
		for (int p = 0; p < cfg.swarmsize_pp; p++)
		{
			if (Fposition[p] < Fpbest[p])
			{
				Fpbest[p] = Fposition[p];
				for (int d = 0; d < cfg.xdim; d++)
				{
					pbest[d][p] = position[d][p];
				}
			}
		}
		// Find index of lowest pbest. This sadly has to be done the C way because I can't into C++
		// Also normal arrays are more sympathetic
		index_argmin = 0;
		double min = Fpbest[0];
		for (int p = 1; p < cfg.swarmsize_pp; p++)
		{
			if (Fpbest[p] < min)
			{
				min = Fpbest[p];
				index_argmin = p;
			}
		}
		// Assign gbest from best pbest
		Fgbest = Fpbest[index_argmin];
		for (int d = 0; d < cfg.xdim; d++)
		{
			gbest[d] = pbest[d][index_argmin];
		}
		// Calculate velocities and new positions combined with random variation
		for (int d = 0; d < cfg.xdim; d++)
		{
			for (int p = 1; p < cfg.swarmsize_pp; p++)
			{
				position[d][p] = 	position[d][p]
									+ cfg.c_inertia * velocity[d][p]
									+ cfg.c_personal * variation(gen) * (pbest[d][p] - position[d][p])
									+ cfg.c_group * variation(gen) * (gbest[d] - position[d][p]);
			}
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
	double Fposition[cfg.swarmsize_bif] = {0};
	double Fpbest[cfg.swarmsize_bif];
	for (int p = 0; p < cfg.swarmsize_bif; p++)
	{
		Fpbest[p] = std::numeric_limits<double>::max();
	}
	double Fgbest = std::numeric_limits<double>::max();
	double *gbest = new double[cfg.ldim]; // Maybe a shared pointer in the future?
	double position[cfg.ldim][cfg.swarmsize_bif];
	double pbest[cfg.ldim][cfg.swarmsize_bif];
	double velocity[cfg.ldim][cfg.swarmsize_bif] = {0};
	// Assign random position one by one within bounds. May need optimization.
	for (int d = 0; d < cfg.xdim; d++)
	{
		std::uniform_real_distribution<> searchspace(cfg.xmin[d], cfg.xmax[d]);
		for (int p = 0; p < cfg.swarmsize_bif; p++)
		{
			position[d][p] = searchspace(gen);
		}
		std::copy(&position[d][0], &position[d][cfg.swarmsize_bif], pbest[0]); // Save memcpy z to p
	}
	// End of Initialization

	std::uniform_real_distribution<> variation(0.0, 1.0); // Used later in next position calculation
	while ((iterations < cfg.iterations_bif) && (Fgbest > cfg.Cstop_bif))
	{
		iterations++;
		bool isOutOfBound = false;
		// Check for out of bound positions and don't apply Fpp to them
		for (int p = 0; p < cfg.swarmsize_bif; p++)
		{
			for (int d = 0; d < cfg.ldim; d++)
			{
				if ((position[d][p] < cfg.lmin[d]) || (position[d][p] > cfg.lmax[d]))
				{
					isOutOfBound = true;
					break;
				}
			}
			if (!isOutOfBound)
			{
				isOutOfBound = false;
				extrenum_t xp = PSOpp(cfg, position[p]);
				Fposition[p] = cfg.objective_bif(cfg, xp, position[p]);
				delete[] xp.point; // This memory gets called in PSOpp!
			}
		}
		// Check if value of position is lower than pbest best and overwrite them if yes
		for (int p = 0; p < cfg.swarmsize_bif; p++)
		{
			if (Fposition[p] < Fpbest[p])
			{
				Fpbest[p] = Fposition[p];
				for (int d = 0; d < cfg.ldim; d++)
				{
					pbest[d][p] = position[d][p];
				}
			}
		}
		// Find index of lowest pbest
		index_argmin = 0;
		double min = Fpbest[0];
		for (int p = 1; p < cfg.swarmsize_bif; p++)
		{
			if (Fpbest[p] < min)
			{
				min = Fpbest[p];
				index_argmin = p;
			}
		}
		// Assign gbest from best pbest
		Fgbest = Fpbest[index_argmin];
		for (int d = 0; d < cfg.ldim; d++)
		{
			gbest[d] = pbest[d][index_argmin];
		}
		// Calculate velocities and new positions combined with random variation
		for (int d = 0; d < cfg.ldim; d++)
		{
			for (int p = 1; p < cfg.swarmsize_bif; p++)
			{
				position[d][p] = 	position[d][p]
									+ cfg.c_inertia * velocity[d][p]
									+ cfg.c_personal * variation(gen) * (pbest[d][p] - position[d][p])
									+ cfg.c_group * variation(gen) * (gbest[d] - position[d][p]);
			}
		}
	}

	extrenum_t best
	{
		.value = Fgbest,
		.point = gbest
	};

	return best; // Don't forget delete[] outside this call! (:
}