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

#if defined(_WIN32)
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <stdlib.h>

#include "config.h"
#include "err.h"
#include "init.h"
#include "potential.h"

/*
 * Allocate position array and generate positions according to a
 * simple lattice.
 */
int init_positions(double ***pos, struct config *conf)
{
    int i, j, k, n, atoms, atoms_per_dimension;
    double box, grid_size;
    double *p;

    atoms = conf->atoms;
    *pos = malloc(sizeof(**pos) * atoms);
    if (*pos) {
        for (i = 0; i < atoms; ++i) {
            (*pos)[i] = malloc(sizeof(***pos) * DIM);
            if (!(*pos)[i]) {
                free(*pos);
                return E_MALLOC;
            }
        }
    } else {
        return E_MALLOC;
    }

    box = conf->box;
    atoms_per_dimension = (int)(pow((double)atoms, (1.0 / DIM)) + 1.);
    grid_size = box / (double)atoms_per_dimension;

    n = 0;
    for (i = 0; i < atoms_per_dimension; ++i) {
        for (j = 0; j < atoms_per_dimension; ++j) {
            for (k = 0; k < atoms_per_dimension; ++k) {
                if (n >= atoms)
                    return E_SUCCESS;
                p = (*pos)[n];
                p[X] = i * grid_size;
                p[Y] = j * grid_size;
                p[Z] = k * grid_size;
                ++n;
            }
        }
    }

    return E_GENERAL;
}

/*
 * TODO: Set the velocities to have a magnitude consistent with a
 * temperature stored in the config struct.
 */
/*
 * Allocate and assign initial velocities using random numbers.  The
 * magnitude should be distributed according to a Maxwell-Boltzmann
 * distribution corresponding to a desired temperature.  Velocities
 * are adjusted according to the center of momentum.
 */
int init_velocities(double ***vel, double *mass, struct config *conf)
{
    int i, atoms;
    double *v;

    atoms = conf->atoms;
    *vel = malloc(sizeof(**vel) * atoms);
    if (*vel) {
        for (i = 0; i < atoms; ++i) {
            (*vel)[i] = malloc(sizeof(***vel) * DIM);
            if (!(*vel)[i]) {
                free(*vel);
                return E_MALLOC;
            }
        }
    } else {
        return E_MALLOC;
    }

    for (i = 0; i < atoms; ++i) {
        v = (*vel)[i];
        v[X] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        v[Y] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        v[Z] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
    }

    scale_velocities(*vel, mass, conf);

    return E_SUCCESS;
}

void scale_velocities(double **vel, double *mass, struct config *conf)
{
    int i;
    double atoms, ke, temperature, scale;
    double *v;
    double centroid[DIM] = { 0.0, 0.0, 0.0 };

    atoms = conf->atoms;
    for (i = 0;    i < atoms; ++i) {
        v = vel[i];
        centroid[X] += mass[i] * v[X];
        centroid[Y] += mass[i] * v[Y];
        centroid[Z] += mass[i] * v[Z];
    }

    centroid[X] /= atoms;
    centroid[Y] /= atoms;
    centroid[Z] /= atoms;

    ke = 0.0;
    for (i = 0; i < atoms; ++i) {
        v = vel[i];
        v[X] -= centroid[X];
        v[Y] -= centroid[Y];
        v[Z] -= centroid[Z];
        ke += mass[i] * (v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
    }

    ke *= 0.5;
    temperature = conf->temperature * KB;
    scale = sqrt(1.5 * atoms * temperature / ke);

    for (i = 0; i < atoms; ++i) {
        v = vel[i];
        v[X] *= scale;
        v[Y] *= scale;
        v[Z] *= scale;
    }
}

/*
 * Allocate and initialize the accerelation array by calculating the
 * forces betwen the atoms at their initial positions.
 */
int init_accelerations(double ***acc, double **pos, double *mass,
                       struct config *conf)
{
    int i, atoms;
    double *a;

    atoms = conf->atoms;
    *acc = malloc(sizeof(**acc) * atoms);
    if (*acc) {
        for (i = 0; i < atoms; ++i) {
            (*acc)[i] = malloc(sizeof(***acc) * DIM);
            if (!(*acc)[i]) {
                free(*acc);
                return E_MALLOC;
            }
        }
    } else {
        return E_MALLOC;
    }

    /* Put the forces into acc. */
    forces(pos, conf, *acc);

    for (i = 0; i < atoms; ++i) {
        a = (*acc)[i];
        a[X] /= mass[i];
        a[Y] /= mass[i];
        a[Z] /= mass[i];
    }

    return E_SUCCESS;
}

/*
 * Properly free a two-dimensional vector with first dimension of
 * length `size'.
 */
void destroy_vector(double **vec, int size)
{
    int i;

    for (i = 0; i < size; ++i)
        free(vec[i]);
    free(vec);
}
