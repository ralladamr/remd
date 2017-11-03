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

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "config.h"

void md(double **pos, double **vel, double **acc, double *mass,
        struct config *conf);
void velocity_verlet(double **pos, double **vel, double **acc, double *mass,
                     struct config *conf);
double virial(double **pos, double **acc, double *mass, double temp,
              struct config *conf);

#endif
