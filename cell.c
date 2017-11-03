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

#include <stdlib.h>
#include <math.h>

#include "cell.h"
#include "config.h"
#include "err.h"

/* The length of each cell in the list. */
static double cell_length;

/* The number of cells per dimension in the list. */
static int cells_per_dim;

/* The total number of cells in the list. */
static int total_cells;

/* The list of neighboring cells. */
static int **neighbors;

static void check_neighbors(int *n);

/*
 * Initialize the cell list. For now, we assume each cell has length
 * equal to the cutoff radius of the simulation.
 */
int init_cell_list(struct cell ***cell_list, double **pos, struct config *conf)
{
    double box, cut;
    int i;

    box = conf->box;
    cut = conf->cutoff;
    cell_length = box / (int)(box / cut);
    cells_per_dim = (int)(box / cell_length);
    total_cells = cells_per_dim;
    for (i = 0; i < DIM - 1; ++i) {
        total_cells *= cells_per_dim;
    }

    *cell_list = malloc(sizeof(*cell_list) * total_cells);
    if (*cell_list) {
        for (i = 0; i < total_cells; ++i) {
            (*cell_list)[i] = malloc(sizeof(**cell_list));
            (*cell_list)[i]->head = NULL;
        }
    } else {
        return E_MALLOC;
    }

    update_cell_list(*cell_list, pos, conf->atoms);

    return E_SUCCESS;
}

void update_cell_list(struct cell **cell_list, double **pos, int atoms)
{
    int i, ix, iy, iz, index;
    double *p;
    struct cell *cur;
    struct node *new;

    for (i = 0; i < atoms; ++i) {
        p = pos[i];
        ix = (int)(p[X] / cell_length);
        iy = (int)(p[Y] / cell_length);
        iz = (int)(p[Z] / cell_length);
        index = ix + iy * cells_per_dim + iz * cells_per_dim * cells_per_dim;
        cur = cell_list[index];
        new = malloc(sizeof(*new));
        new->atom = i;
        new->next = cur->head;
        cur->head = new;
    }
}

/* Allocate and initialize an array of cell neighbors base on a cell list.*/
int init_neighbors(struct cell **cell_list)
{
    int i;
    int *n;

    neighbors = malloc(sizeof(neighbors) * total_cells);
    if (*neighbors) {
        for (i = 0; i < total_cells; ++i) {
            neighbors[i] = malloc(sizeof(*neighbors) * FACES);
            if (!neighbors[i]) {
                free(*neighbors);
                return E_MALLOC;
            }
        }
    }
    else {
        return E_MALLOC;
    }

    for (i = 0; i < total_cells; ++i) {
        n = neighbors[i];
        n[0] = i + 1;
        n[1] = i - 1;
        n[2] = i + cells_per_dim;
        n[3] = i - cells_per_dim;
        n[4] = i + cells_per_dim + 1;
        n[5] = i + cells_per_dim - 1;
        n[6] = i - cells_per_dim + 1;
        n[7] = i - cells_per_dim - 1;
        n[8] = i + cells_per_dim * cells_per_dim;
        n[9] = i + cells_per_dim * cells_per_dim + 1;
        n[10] = i + cells_per_dim * cells_per_dim - 1;
        n[11] = i + cells_per_dim * (cells_per_dim + 1);
        n[11] = i + cells_per_dim * (cells_per_dim - 1);
        n[12] = i + cells_per_dim * (cells_per_dim + 1) + 1;
        n[13] = i + cells_per_dim * (cells_per_dim + 1) - 1;
        n[15] = i + cells_per_dim * (cells_per_dim - 1) + 1;
        n[16] = i + cells_per_dim * (cells_per_dim - 1) - 1;
        n[17] = i - cells_per_dim * cells_per_dim;
        n[18] = i - cells_per_dim * cells_per_dim + 1;
        n[19] = i - cells_per_dim * cells_per_dim - 1;
        n[20] = i - cells_per_dim * (cells_per_dim + 1);
        n[21] = i - cells_per_dim * (cells_per_dim - 1);
        n[22] = i - cells_per_dim * (cells_per_dim + 1) + 1;
        n[23] = i - cells_per_dim * (cells_per_dim + 1) - 1;
        n[24] = i - cells_per_dim * (cells_per_dim - 1) + 1;
        n[25] = i - cells_per_dim * (cells_per_dim - 1) - 1;
        check_neighbors(n);
    }

    return E_SUCCESS;
}

static void check_neighbors(int n[FACES])
{
    int i, ix, iy, iz;

    for (i = 0; i < FACES; ++i) {
        iz = n[i] / (cells_per_dim * cells_per_dim);
        iy = n[i] - (cells_per_dim * cells_per_dim) / cells_per_dim;
        ix = n[i] - (cells_per_dim * (cells_per_dim - iy))
            / cells_per_dim;

        if (ix < 0) {
            ix = cells_per_dim;
        }

        if (ix > cells_per_dim) {
            ix = 0;
        }

        if (iy < 0) {
            iy = cells_per_dim;
        }

        if (iy > cells_per_dim) {
            iy = 0;
        }

        if (iz < 0) {
            iz = cells_per_dim;
        }

        if (iz > cells_per_dim) {
            iz = 0;
        }

        n[i] = ix + iy * cells_per_dim
            + iz * cells_per_dim * cells_per_dim;
    }
}

void get_neighbors(int i, int n[FACES])
{
    n = neighbors[i];
}


int get_total_cells()
{
    return total_cells;
}

void destroy_cell_list(struct cell **cell_list)
{
    int i;
    struct cell *cur;

    for (i = 0; i < total_cells; ++i) {
        cur = cell_list[i];
        if (cur->head) {
            free(cur->head);
        }
        free(cur);
    }
    free(cell_list);
}

void destroy_neighbor_list()
{
    int i;

    for (i = 0; i < total_cells; ++i) {
        free(neighbors[i]);
    }
    free(neighbors);
}
