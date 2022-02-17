// © 2022 Robert Hoffmann <robert.hoffmann@smail.emt.h-brs.de>
// I'll release this under a license once I decided which.
#include "nlpso.hpp"

double* PSOpp(nlpso_cfg_t cfg, double *lambda)
{
	// Seeding rand function. I just copied this from cppreference lol
	std::random_device rd;
	std::mt19937 gen(rd());
	// Start of Initialization
	int index_argmin;
	int iterations = 0;
	double Fposition[cfg.swarmsize_pp] = {0};
	double Fpbest[cfg.swarmsize_pp] = {std::numeric_limits<double>::max()};
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
			std::cout << position[d][p] << "  ";// Remove after debugging //
		}
		std::copy(&position[d][0], &position[d][16], pbest[0]); // Save memcpy z to p
		std::cout << std::endl;// Remove after debugging //
	}
	// End of Initialization

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
				Fposition[p] = cfg.weight_pp(position[0][p], lambda, cfg.period);
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
		// Find lowest pbest best. This sadly has to be done the C way because I can't into C++
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
			gbest[d] = pbest[d][index_argmin];// May not needed!
		}
		// Calculate velocities and new positions combined with random variation
		std::uniform_real_distribution<> variation(0.0, 1.0);
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

	return gbest; // Don't forget delete[] outside this call! (:
}