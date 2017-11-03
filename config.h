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

#ifndef CONFIG_H
#define CONFIG_H

#define DIM 3

/* Meters in one angstrom */
#define ANGSTROM 1.0e-10

/* Boltzman constant */
#define KB 1.3806488e-23

/* Kg in one amu */
#define AMU 1.660468e-27

#define X 0
#define Y 1
#define Z 2

struct config {
    int atoms;
    long steps;
    int lorentz_berthelot;
    double time_step;
    double cutoff;
    double box;
    double temperature;
    int write_steps;
    int equilibration_steps;
};

struct config read_config(char *config_file);

#endif