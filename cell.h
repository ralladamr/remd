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

#ifndef CELL_H
#define CELL_H

#include "config.h"

#define FACES 26

struct node {
    int atom;
    struct node *next;
};

struct cell {
    struct node *head;
};

int init_cell_list(struct cell ***cell_list, double **pos, struct config *conf);
void update_cell_list(struct cell **cell_list, double **pos, int atoms);
int init_neighbors(struct cell **cell_list);
void get_neighbors(int i, int n[FACES]);
int get_total_cells();
void destroy_cell_list(struct cell **cell_list);

#endif

