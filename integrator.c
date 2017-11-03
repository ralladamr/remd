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

#include "integrator.h"
#include "potential.h"

/*
 * The main molecular dynamics loop.
 */
void md(double **pos, double **vel, double **acc, double *mass,
        struct config *conf)
{
    int atoms, steps, equilibration, write;
    double avg_ke, ke, avg_pe, pe, avg_temp, temp;
    double avg_pressure, pressure, initial_energy;
    long i;

    atoms = conf->atoms;
    ke = kinetic_energy(vel, mass, atoms);
    pe = potential_energy(pos, conf);
    initial_energy = ke + pe;

    avg_ke = 0.0;
    avg_pe = 0.0;
    avg_temp = 0.0;
    avg_pressure = 0.0;
    steps = conf->steps;
    write = conf->write_steps;
    equilibration = conf->equilibration_steps;
    for (i = 0; i < steps; ++i) {
        velocity_verlet(pos, vel, acc, mass, conf);

        if ((i + 1) % write == 0 && (i + 1) > equilibration) {
            ke = kinetic_energy(vel, mass, atoms);
            avg_ke += ke;
            temp = 2.0 * (ke / KB) / 3.0 / atoms;
            avg_temp += temp;
            pe = potential_energy(pos, conf);
            avg_pe += pe;
            pressure = virial(pos, acc, mass, temp, conf);
            avg_pressure += pressure;

            printf("%e %e %e %f %f %e\n", ke, pe, ke + pe,
                   (ke + pe - initial_energy) / initial_energy, temp,
                   pressure);
        }
    }

    avg_ke /= (steps - equilibration) / write;
    avg_pe /= (steps - equilibration) / write;
    avg_temp /= (steps - equilibration) / write;
    avg_pressure /= (steps - equilibration) / write;

    printf("\nAVERAGES\n");
    printf("========\n");
    printf("%e %e %f %e\n", avg_ke, avg_pe, avg_temp, avg_pressure);
}

/*
 * Perform velocity Verlet integration using the positions,
 * velocities, accelerations and masses passed as pos, vel, acc, and
 * mass, respectively.
 */
void velocity_verlet(double **pos, double **vel, double **acc, double *mass,
                     struct config *conf)
{
    int i, atoms;
    double dt, box;
    double *p, *v, *a;

    dt = conf->time_step;
    box = conf->box;
    atoms = conf->atoms;
#pragma omp parallel for shared(atoms, pos, vel, acc, dt, box) \
    private(i, p, v, a)
    for (i = 0; i < atoms; ++i) {
        p = pos[i];
        v = vel[i];
        a = acc[i];
        p[X] += dt * v[X] + 0.5 * dt * dt * a[X];
        p[X] -= floor(p[X] / box) * box;
        p[Y] += dt * v[Y] + 0.5 * dt * dt * a[Y];
        p[Y] -= floor(p[Y] / box) * box;
        p[Z] += dt * v[Z] + 0.5 * dt * dt * a[Z];
        p[Z] -= floor(p[Z] / box) * box;
        v[X] += 0.5 * dt * a[X];
        v[Y] += 0.5 * dt * a[Y];
        v[Z] += 0.5 * dt * a[Z];
    }

    forces(pos, conf, acc);

#pragma omp parallel for shared(atoms, acc, mass) private(i, a)
    for (i = 0; i < atoms; ++i) {
        a = acc[i];
        a[X] /= mass[i];
        a[Y] /= mass[i];
        a[Z] /= mass[i];
    }

#pragma omp parallel for shared(atoms, vel, acc, dt) private(i, a, v)
    for (i = 0; i < atoms; ++i) {
        v = vel[i];
        a = acc[i];
        v[X] += 0.5 * dt * a[X];
        v[Y] += 0.5 * dt * a[Y];
        v[Z] += 0.5 * dt * a[Z];
    }
}

/*
 * Calculate the pressure using the virial function.
 */
double virial(double **pos, double **acc, double *mass, double temp,
              struct config *conf)
{
    int atoms, i;
    double box, volume, m, pressure;
    double *a, *p;

    atoms = conf->atoms;
    box = conf->box;
    volume = box * box * box;
    pressure = 0.0;

#pragma omp parallel for shared(atoms, acc, pos, mass) private(i, a, p, m) \
    reduction(+:pressure)
    for (i = 0; i < atoms; ++i) {
        a = acc[i];
        p = pos[i];
        m = mass[i];
        pressure += p[X] * a[X] / m / DIM;
        pressure += p[Y] * a[Y] / m / DIM;
        pressure += p[Z] * a[Z] / m / DIM;
    }
    pressure /= atoms / volume;
    pressure += atoms * temp * KB / volume;
    return pressure;
}