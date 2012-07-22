#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <math.h>
#include "interp2d.h"

const static double bicubic_inversion_matrix[] = {
    // columns act on:
    // values      x derivs     y derivs      xy derivs
     1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0, // computes a_00
     0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0, // computes a_01
    -3, 3, 0, 0,   0, 0, 0, 0,  -2,-1, 0, 0,   0, 0, 0, 0, // computes a_02
     2,-2, 0, 0,   0, 0, 0, 0,   1, 1, 0, 0,   0, 0, 0, 0, // computes a_03
     
     0, 0, 0, 0,   1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0, // computes a_10
     0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 0, 0, // computes a_11
     0, 0, 0, 0,  -3, 3, 0, 0,   0, 0, 0, 0,  -2,-1, 0, 0, // computes a_12
     0, 0, 0, 0,   2,-2, 0, 0,   0, 0, 0, 0,   1, 1, 0, 0, // computes a_13
     
    -3, 0, 3, 0,  -2, 0,-1, 0,   0, 0, 0, 0,   0, 0, 0, 0, // computes a_20
     0, 0, 0, 0,   0, 0, 0, 0,  -3, 0, 3, 0,  -2, 0,-1, 0, // computes a_21
     9,-9,-9, 9,   6,-6, 3,-3,   6, 3,-6,-3,   4, 2, 2, 1, // computes a_22
    -6, 6, 6,-6,  -4, 4,-2, 2,  -3,-3, 3, 3,  -2,-2,-1,-1, // computes a_23
    
     2, 0,-2, 0,   1, 0, 1, 0,   0, 0, 0, 0,   0, 0, 0, 0, // computes a_30
     0, 0, 0, 0,   0, 0, 0, 0,   2, 0,-2, 0,   1, 0, 1, 0, // computes a_31
    -6, 6, 6,-6,  -3, 3,-3, 3,  -4,-2, 4, 2,  -2,-1,-2,-1, // computes a_32
     4,-4,-4, 4,   2,-2, 2,-2,   2, 2,-2,-2,   1, 1, 1, 1  // computes a_33
};

static int bicubic_init(void* state, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize) {
    // TODO: fill in later if necessary
    return GSL_SUCCESS;
}

static int bicubic_eval(const void* state, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    double xmin, xmax, ymin, ymax;
    // z_values_and_derivatives has the following structure:
    //  [0] = value of the function at (xmin, ymin)
    //  [1] = value of the function at (xmin, ymax)
    //  [2] = value of the function at (xmax, ymin)
    //  [3] = value of the function at (xmax, ymax)
    //  [4] = x derivative of the function at (xmin, ymin)
    //  [5] = x derivative of the function at (xmin, ymax)
    //  [6] = x derivative of the function at (xmax, ymin)
    //  [7] = x derivative of the function at (xmax, ymax)
    //  [8] = y derivative of the function at (xmin, ymin)
    //  ...
    //  [12] = xy derivative of the function at (xmin,ymin)
    //  ...
    double z_values_and_derivatives[16];
    // coefficients has the following structure:
    //  [0] = a_00
    //  [1] = a_01
    //  [2] = a_02
    //  [3] = a_03
    //  [4] = a_10
    //  ...
    // The final value comes from summing a_ij t^i u^j (see t and u below)
    double coefficients[16];
    // dx and dy are the size of the grid cell
    double dx, dy;
    // t and u are the positions within the grid cell at which we are computing
    // the interpolation, in units of grid cell size
    double t, u;
    double z_accumulator;
    size_t xi, yi;
    size_t i, j;
    // First compute the indices into the data arrays where we are interpolating
    if (xa != NULL) {
        xi = gsl_interp_accel_find(xa, xarr, xsize, x);
    }
    else {
        xi = gsl_interp_bsearch(xarr, x, 0, xsize - 1);
    }
    if (ya != NULL) {
        yi = gsl_interp_accel_find(ya, yarr, ysize, y);
    }
    else {
        yi = gsl_interp_bsearch(yarr, y, 0, ysize - 1);
    }
    // Find the minimum and maximum values on the grid cell in each dimension
    xmin = xarr[xi];
    xmax = xarr[xi + 1];
    ymin = yarr[yi];
    ymax = yarr[yi + 1];
    z_values_and_derivatives[0] = zarr[INDEX_2D(xi, yi, xsize, ysize)];
    z_values_and_derivatives[1] = zarr[INDEX_2D(xi, yi + 1, xsize, ysize)];
    z_values_and_derivatives[2] = zarr[INDEX_2D(xi + 1, yi, xsize, ysize)];
    z_values_and_derivatives[3] = zarr[INDEX_2D(xi + 1, yi + 1, xsize, ysize)];
    // Get the width and height of the grid cell
    dx = xmax - xmin;
    dy = ymax - ymin;
    // Determine the values of the function's derivatives at each grid point
    // using a 3-point numerical differentiation scheme
    z_values_and_derivatives[4]  = 0.5 / dx * (
          zarr[INDEX_2D(xi + 1, yi, xsize, ysize)]
        - zarr[INDEX_2D(xi - 1, yi, xsize, ysize)]);
    z_values_and_derivatives[5]  = 0.5 / dx * (
          zarr[INDEX_2D(xi + 1, yi + 1, xsize, ysize)]
        - zarr[INDEX_2D(xi - 1, yi + 1, xsize, ysize)]);
    z_values_and_derivatives[6]  = 0.5 / dx * (
          zarr[INDEX_2D(xi + 2, yi, xsize, ysize)]
        - zarr[INDEX_2D(xi, yi, xsize, ysize)]);
    z_values_and_derivatives[7]  = 0.5 / dx * (
          zarr[INDEX_2D(xi + 2, yi + 1, xsize, ysize)]
        - zarr[INDEX_2D(xi, yi + 1, xsize, ysize)]);
    z_values_and_derivatives[8]  = 0.5 / dy * (
          zarr[INDEX_2D(xi, yi + 1, xsize, ysize)]
        - zarr[INDEX_2D(xi, yi - 1, xsize, ysize)]);
    z_values_and_derivatives[9]  = 0.5 / dy * (
          zarr[INDEX_2D(xi, yi + 2, xsize, ysize)]
        - zarr[INDEX_2D(xi, yi, xsize, ysize)]);
    z_values_and_derivatives[10] = 0.5 / dy * (
          zarr[INDEX_2D(xi + 1, yi + 1, xsize, ysize)]
        - zarr[INDEX_2D(xi + 1, yi - 1, xsize, ysize)]);
    z_values_and_derivatives[11] = 0.5 / dy * (
          zarr[INDEX_2D(xi + 1, yi + 2, xsize, ysize)]
        - zarr[INDEX_2D(xi + 1, yi, xsize, ysize)]);
    z_values_and_derivatives[12] = 0.25 / (dx * dy) * (
          zarr[INDEX_2D(xi + 1, yi + 1, xsize, ysize)]
        - zarr[INDEX_2D(xi + 1, yi - 1, xsize, ysize)]
        - zarr[INDEX_2D(xi - 1, yi + 1, xsize, ysize)]
        + zarr[INDEX_2D(xi - 1, yi - 1, xsize, ysize)]);
    z_values_and_derivatives[13] = 0.25 / (dx * dy) * (
          zarr[INDEX_2D(xi + 1, yi + 2, xsize, ysize)]
        - zarr[INDEX_2D(xi + 1, yi,     xsize, ysize)]
        - zarr[INDEX_2D(xi - 1, yi + 2, xsize, ysize)]
        + zarr[INDEX_2D(xi - 1, yi,     xsize, ysize)]);
    z_values_and_derivatives[14] = 0.25 / (dx * dy) * (
          zarr[INDEX_2D(xi + 2, yi + 1, xsize, ysize)]
        - zarr[INDEX_2D(xi + 2, yi - 1, xsize, ysize)]
        - zarr[INDEX_2D(xi,     yi + 1, xsize, ysize)]
        + zarr[INDEX_2D(xi,     yi - 1, xsize, ysize)]);
    z_values_and_derivatives[15] = 0.25 / (dx * dy) * (
          zarr[INDEX_2D(xi + 2, yi + 2, xsize, ysize)]
        - zarr[INDEX_2D(xi + 2, yi,     xsize, ysize)]
        - zarr[INDEX_2D(xi,     yi + 2, xsize, ysize)]
        + zarr[INDEX_2D(xi,     yi,     xsize, ysize)]);
    // and the relative position within the cell where we are interpolating
    t = (x - xmin)/dx;
    u = (y - ymin)/dy;
    // Matrix multiplication of the z values and derivatives with the 
    // constant matrix above to produce the interpolation coefficients
    for (i = 0; i < 16; i++) {
        coefficients[i] = 0;
        for (j = 0; j < 16; j++) {
            coefficients[i] += bicubic_inversion_matrix[16*i+j] * z_values_and_derivatives[j];
        }
    }
    // Compute the new value as sum(a_ij t^i u^j)
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            z_accumulator += coefficients[4*i+j] * pow(t, i) * pow(u, j);
        }
    }
    *z = z_accumulator;
    return GSL_SUCCESS;
}

static const interp2d_type bicubic_type = {
    "bicubic",
    4,
    NULL,
    &bicubic_init,
    &bicubic_eval
};

const interp2d_type* interp2d_bicubic = &bicubic_type;
