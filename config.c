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

#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "err.h"

static void check_read(int err);

/* Read a set of system parameters from a config file. */
struct config read_config(char *config_file)
{
    struct config conf;
    FILE *fp;
    int err;

#if defined(_WIN32)
    fopen_s(&fp, config_file, "r");
#else
    fp = fopen(config_file, "r");
#endif

    if (fp) {
#if defined(_WIN32)
        err = fscanf_s(fp, "%d", &conf.atoms);
        check_read(err);
        err = fscanf_s(fp, "%ld", &conf.steps);
        check_read(err);
        err = fscanf_s(fp, "%d", &conf.equilibration_steps);
        check_read(err);
        err = fscanf_s(fp, "%d", &conf.write_steps);
        check_read(err);
        err = fscanf_s(fp, "%lf", &conf.time_step);
        check_read(err);
        err = fscanf_s(fp, "%lf", &conf.box);
        check_read(err);
        err = fscanf_s(fp, "%lf", &conf.cutoff);
        check_read(err);
        err = fscanf_s(fp, "%lf", &conf.temperature);
#else
        err = fscanf(fp, "%d", &conf.atoms);
        check_read(err);
        err = fscanf(fp, "%ld", &conf.steps);
        check_read(err);
        err = fscanf(fp, "%d", &conf.equilibration_steps);
        check_read(err);
        err = fscanf(fp, "%d", &conf.write_steps);
        check_read(err);
        err = fscanf(fp, "%lf", &conf.time_step);
        check_read(err);
        err = fscanf(fp, "%lf", &conf.box);
        check_read(err);
        err = fscanf(fp, "%lf", &conf.cutoff);
        check_read(err);
        err = fscanf(fp, "%lf", &conf.temperature);
        check_read(err);
#endif
        fclose(fp);
    } else {
        fprintf(stderr, "Configuration file not found!\n");
        exit(EXIT_FAILURE);
    }

    conf.lorentz_berthelot = 1;

    return conf;
}

static void check_read(int err)
{
    if (EOF == err || 0 == err) {
        fprintf(stderr, "Could not read all configuration parameters!\n");
        exit(EXIT_FAILURE);
    }
}
