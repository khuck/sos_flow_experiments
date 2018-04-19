/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include "globals.h"
#include "simple_timer.h"
#include "mpi.h"
#include "adios.h"

/*************************************************************/
/*          Example of writing matrix in ADIOS               */
/*                                                           */
/*************************************************************/
void do_adios(Matrix<double> &matrix)
{
    char        filename [256];
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;
    int rows = matrix.Rows();
    int cols = matrix.Columns();

    strcpy (filename, "matrix.bp");
    simple_timer t("ADIOS");
    adios_open (&adios_handle, "matrix", filename, "w", MPI_COMM_WORLD);
    adios_groupsize = 4 + 4 + 8 * (rows) * (cols);
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "NX", &rows);
    adios_write (adios_handle, "NY", &cols);
    adios_write (adios_handle, "var_double_2Darray", matrix.Data());
    adios_close (adios_handle);
    MPI_Barrier (MPI_COMM_WORLD);
    return;
}
