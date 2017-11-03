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
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "err.h"
#include "potential.h"

#define NINT(x) ((x) < 0 ? ((long)(x - 0.5)) : ((long)(x + 0.5)))

/*
 * This array holds the N x N array of Lennard-Jones interaction
 * energies.
 */
static double **eps;

/* This array holds the N x N array of Lennard-Jones distances. */
static double **sig;

/* Free the sig and eps arrays properly. */
void destroy_potential_params(int size)
{
    int i;
    for (i = 0; i < size; ++i) {
        free(sig[i]);
        free(eps[i]);
    }
    free(sig);
    free(eps);
}

/*
 * Malloc and initialize the sig and eps arrays.  The sigma and
 * epsilon pointers passed as parameters contain the like-like
 * interactions; the Lorentz-Berthelot mixing rules are used to
 * generate the unlike pairs.
 */
int init_potential_params(double *sigma, double *epsilon, struct config *conf)
{
    int  atoms, i, j;

    atoms = conf->atoms;
    sig = malloc(sizeof(*sig) * atoms);
    if (sig) {
        for (i = 0; i < atoms; ++i) {
            sig[i] = malloc(sizeof *sig[i] * atoms);
            if (!sig[i]) {
                return E_MALLOC;
            }
        }
    } else {
        return E_MALLOC;
    }

    eps = malloc(sizeof(*eps) * atoms);
    if (eps) {
        for (i = 0; i < atoms; ++i) {
            eps[i] = malloc(sizeof *eps[i] * atoms);
            if (!eps[i]) {
                return E_MALLOC;
            }
        }
    } else {
        return E_MALLOC;
    }

    for (i = 0; i < atoms; ++i) {
        sig[i][i] = sigma[i];
        eps[i][i] = epsilon[i];
    }

    /* Only LB is supported for now! */
    if (conf->lorentz_berthelot) {
        for (i = 0; i < atoms; ++i) {
            for (j = 0; j < atoms; ++j) {
                if (i == j) {
                    continue;
                }
                sig[i][j] = (sigma[i] + sigma[j]) / 2;
                eps[i][j] = sqrt(epsilon[i] * epsilon[j]);
            }
        }
    }

    return E_SUCCESS;
}

/*
 * Compute the forces between each atom in the system.  Place the
 * results in the frc array.
 */
void forces(double **pos, struct config *conf, double **frc)
{
    int i, j, atoms;
    double box, cut, dx, dy, dz, r2, r6, rc2, rc6, sigij2, fij;
    double *p1, *p2, *f1, *f2;

    atoms = conf->atoms;
    for (i = 0; i < atoms; ++i) {
        f1 = frc[i];
        f1[X] = 0.0;
        f1[Y] = 0.0;
        f1[Z] = 0.0;
    }

    box = conf->box;
    cut = conf->cutoff;
    rc2 = cut * cut;

#pragma omp parallel for shared(frc, pos, box, rc2, sig, eps, atoms) \
    private(i, j, p1, p2, f1, f2, dx, dy, dz, r2, r6, rc6, sigij2, fij)
    for (i = 0; i < atoms; ++i) {
        p1 = pos[i];
        f1 = frc[i];
        for (j = i + 1; j < atoms; ++j) {
            p2 = pos[j];
            f2 = frc[j];
            dx = (p1[X] - p2[X]);
            dx -= NINT(dx / box) * box;
            dy = (p1[Y] - p2[Y]);
            dy -= NINT(dy / box) * box;
            dz = (p1[Z] - p2[Z]);
            dz -= NINT(dz / box) * box;
            r2 = dx * dx + dy * dy + dz * dz;

            if (r2 < rc2) {
                sigij2 = sig[i][j];
                sigij2 *= sigij2;
                r6 = sigij2 / r2;
                r6 *= r6 * r6;
                rc6 = sigij2 / rc2;
                rc6 *= rc6 * rc6;
                fij = (24.0 * eps[i][j] * (r6 * (2.0 * r6 - 1.0) / r2
                       - rc6 * (2.0 * rc6 - 1.0) / rc2));
                f1[X] += fij * dx;
                f1[Y] += fij * dy;
                f1[Z] += fij * dz;
                f2[X] -= fij * dx;
                f2[Y] -= fij * dy;
                f2[Z] -= fij * dz;
            }
        }
    }
}

/*
 * Compute and return the total Lennard-Jones potential energy of the
 * atoms in the system.
 */
double potential_energy(double **pos, struct config *conf)
{
    int i, j, atoms;
    double box, cut, dx, dy, dz, r2, r6, rc2, rc6, sigij2, energy;
    double *p1, *p2;

    atoms = conf->atoms;
    box = conf->box;
    cut = conf->cutoff;
    rc2 = cut * cut;

    energy = 0.0;

#pragma omp parallel for shared(atoms, pos, rc2, sig, eps) \
    private(i, j, p1, p2, dx, dy, dz, r2, r6, rc6, sigij2) reduction (+:energy)
    for (i = 0; i < atoms - 1; ++i) {
        p1 = pos[i];
        for (j = i + 1; j < atoms; ++j) {
            p2 = pos[j];
            dx = (p1[X] - p2[X]);
            dx -= NINT(dx / box) * box;
            dy = (p1[Y] - p2[Y]);
            dy -= NINT(dy / box) * box;
            dz = (p1[Z] - p2[Z]);
            dz -= NINT(dz / box) * box;
            r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc2) {
                sigij2 = sig[i][j];
                sigij2 *= sigij2;
                r6 = sigij2 / r2;
                r6 *= r6 * r6;
                rc6 = sigij2 / rc2;
                rc6 *= rc6 * rc6;
                energy += (4.0 * eps[i][j] * (r6 * (r6 - 1.0) + 3.0 * rc6
                           * (2.0 * rc6 - 1.0) * r2 / rc2 - rc6
                           * (7.0 * rc6 + 4.0)));
            }
        }
    }
    return energy;
}

/* Compute and return the kinetic energy of the system. */
double kinetic_energy(double **vel, double *mass, int atoms)
{
    int i;
    double energy, m;
    double *v;

    energy = 0.0;

#pragma omp parallel for shared(vel, mass, atoms) private(v, m, i) \
    reduction(+:energy)
    for (i = 0; i < atoms; ++i) {
        v = vel[i];
        m = mass[i];
        energy += m * v[X] * v[X];
        energy += m * v[Y] * v[Y];
        energy += m * v[Z] * v[Z];
    }
    return 0.5 * energy;
}
