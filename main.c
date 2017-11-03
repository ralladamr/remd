/*
 * Copyright (C) 2014 Adam R. Rall
 * This file is part of remd.
 *
 * Remd is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Remd is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with remd; if not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "config.h"
#include "err.h"
#include "init.h"
#include "integrator.h"
#include "potential.h"

/*
 * Remd is a rudimentary example of a molecular dynamics simulation.
 * It performs velocity Verlet integration of a system of
 * Lennard-Jones spheres in a cubic box with periodic boundary
 * conditions.
 */
int main(int argc, char *argv[])
{
    int atoms, i;
    struct config conf;
    double **pos, **vel, **acc;
    double *mass, *sigma, *epsilon;

    /* Begin by initializing the system. */
    srand((unsigned int)time(NULL));
    if (argc > 1) {
        conf = read_config(argv[1]);
    } else {
        fprintf(stderr, "Configuration file required!\n");
        exit(EXIT_FAILURE);
    }
    atoms = conf.atoms;

    /* TODO: These should probably go into separate functions. */
    mass = malloc(sizeof(*mass) * atoms);
    if (!mass) {
        fprintf(stderr, "Malloc failed for mass.\n");
        exit(EXIT_FAILURE);
    }

    sigma = malloc(sizeof(*sigma) * atoms);
    if (!sigma) {
        fprintf(stderr, "Malloc failed for sigma.\n");
        exit(EXIT_FAILURE);
    }

    epsilon = malloc(sizeof(*epsilon) * atoms);
    if (!epsilon) {
        fprintf(stderr, "Malloc failed for epsilon.\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < atoms; ++i) {
        mass[i] = 18 * AMU;
        sigma[i] = 3.40 * ANGSTROM;
        epsilon[i] = 1.6556e-21;
    }

    if (init_potential_params(sigma, epsilon, &conf) != E_SUCCESS) {
        fprintf(stderr, "Malloc failed for potential parameters.\n");
        exit(EXIT_FAILURE);
    }

    free(sigma);
    free(epsilon);

    if (init_positions(&pos, &conf) != E_SUCCESS) {
        fprintf(stderr, "Failed to initialize positions.\n");
        exit(EXIT_FAILURE);
    }

    if (init_velocities(&vel, mass, &conf) != E_SUCCESS) {
        fprintf(stderr, "Failed to initialize velocities.\n");
        exit(EXIT_FAILURE);
    }

    if (init_accelerations(&acc, pos, mass, &conf) != E_SUCCESS) {
        fprintf(stderr, "Failed to initialize accelerations.\n");
        exit(EXIT_FAILURE);
    }

    /* Perform the simulation. */
    md(pos, vel, acc, mass, &conf);

    destroy_vector(pos, atoms);
    destroy_vector(vel, atoms);
    destroy_vector(acc, atoms);
    free(mass);
    destroy_potential_params(atoms);

    return EXIT_SUCCESS;
}
