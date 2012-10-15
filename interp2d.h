/*
 * This file is part of interp2d, a GSL-compatible two-dimensional
 * interpolation library. <http://www.ellipsix.net/interp2d.html>
 * 
 * Copyright 2012 David Zaslavsky
 * Portions based on GNU GSL interpolation code,
 *  copyright 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __INTERP_2D_H__
#define __INTERP_2D_H__

#include <gsl/gsl_interp.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Represents a 2D interpolation algorithm. Example: linear interpolation.
 */
typedef struct {
    /** The name of the algorithm. Use interp2d_name to access this field. */
    const char* name;
    /** The minimum number of points in each dimension required by the algorithm. Use interp2d_min_size to access this field. */
    unsigned int min_size;
    /** The method that allocates memory for an interpolation object of this type. */
    void* (*alloc)(size_t xsize, size_t ysize);
    /** The method that initializes the interpolation type. */ 
    int (*init)(void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize);
    /** The method that evaluates the interpolation at the given point. */
    int (*eval)(const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z);
    /** The method that evaluates the x derivative of the interpolation at the given point. */
    int (*eval_deriv_x) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_p);
    /** The method that evaluates the y derivative of the interpolation at the given point. */
    int (*eval_deriv_y) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_p);
    /** The method that evaluates the second x derivative of the interpolation at the given point. */
    int (*eval_deriv_xx) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_pp);
    /** The method that evaluates the second y derivative of the interpolation at the given point. */
    int (*eval_deriv_xy) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_pp);
    /** The method that evaluates the cross derivative of the interpolation at the given point. */
    int (*eval_deriv_yy) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_pp);
    /** The method that frees the memory. */
    void (*free)(void*);
} interp2d_type;

/**
 * Represents a 2D interpolation instance.
 */
typedef struct {
    /** The type object for the interpolation. */
    const interp2d_type* type;
    /** The minimum value of x for which data have been provided. */
    double xmin;
    /** The maximum value of x for which data have been provided. */
    double xmax;
    /** The minimum value of y for which data have been provided. */
    double ymin;
    /** The maximum value of y for which data have been provided. */
    double ymax;
    /** The number of x values provided. */
    size_t xsize;
    /** The number of y values provided. */
    size_t ysize;
    /** A state object. This is specific to the interpolation type. */
    void* state;
} interp2d;

/** An instance for bilinear interpolation. */
GSL_VAR const interp2d_type* interp2d_bilinear;
// To be added:
// GSL_VAR const interp2d_type* interp2d_bicubic;
// GSL_VAR const interp2d_type* interp2d_bipoly;

/**
 * Return the minimum number of points needed by the given interpolation type.
 */
size_t interp2d_type_min_size(const interp2d_type* T);
/**
 * Return the minimum number of points needed by the type of the given interpolation.
 */
size_t interp2d_min_size(const interp2d* interp);
/**
 * Return the type name of the given interpolation.
 */
const char* interp2d_name(const interp2d* interp);

/**
 * Allocate a new interpolation object. Use this when you want to run an interpolation.
 * 
 * @param T the const struct representing the interpolation algorithm you want to use
 * @param xsize the number of points that will be provided in the x direction
 * @param ysize the number of points that will be provided in the y direction
 * @return the newly allocated object
 */
interp2d* interp2d_alloc(const interp2d_type* T, size_t xsize, size_t ysize);
/**
 * Initialize an interpolation object with data.
 * 
 * This method creates an accelerator cache based on the data provided, but it doesn't actually store the data.
 * 
 * @param interp the interpolation object, previously initialized
 * @param xa the x coordinates of the data, of length xsize
 * @param ya the y coordinates of the data, of length ysize
 * @param za the z coordinates of the data, of length xsize*ysize
 * @param xsize the length of the array xa
 * @param ysize the length of the array ya
 */
int interp2d_init(interp2d* interp, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize);
/**
 * Free the interpolation object.
 */
void interp2d_free(interp2d* interp);

/**
 * Evaluate the interpolating function with the given data.
 * 
 * @param interp the interpolation object, previously initialized
 * @param xarr the x coordinates of the data, of length xsize
 * @param yarr the y coordinates of the data, of length ysize
 * @param zarr the z coordinates of the data, of length xsize*ysize
 * @param xsize the length of the array xa
 * @param ysize the length of the array ya
 * @param xa the accelerator object for the x direction (may be NULL)
 * @param ya the accelerator object for the y direction (may be NULL)
 */
double interp2d_eval(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the x derivative of the interpolating function with the given data.
 * 
 * @param interp the interpolation object, previously initialized
 * @param xarr the x coordinates of the data, of length xsize
 * @param yarr the y coordinates of the data, of length ysize
 * @param zarr the z coordinates of the data, of length xsize*ysize
 * @param xsize the length of the array xa
 * @param ysize the length of the array ya
 * @param xa the accelerator object for the x direction (may be NULL)
 * @param ya the accelerator object for the y direction (may be NULL)
 */
double interp2d_eval_deriv_x(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the y derivative of the interpolating function with the given data.
 * 
 * @param interp the interpolation object, previously initialized
 * @param xarr the x coordinates of the data, of length xsize
 * @param yarr the y coordinates of the data, of length ysize
 * @param zarr the z coordinates of the data, of length xsize*ysize
 * @param xsize the length of the array xa
 * @param ysize the length of the array ya
 * @param xa the accelerator object for the x direction (may be NULL)
 * @param ya the accelerator object for the y direction (may be NULL)
 */
double interp2d_eval_deriv_y(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the second x derivative of the interpolating function with the given data.
 * 
 * @param interp the interpolation object, previously initialized
 * @param xarr the x coordinates of the data, of length xsize
 * @param yarr the y coordinates of the data, of length ysize
 * @param zarr the z coordinates of the data, of length xsize*ysize
 * @param xsize the length of the array xa
 * @param ysize the length of the array ya
 * @param xa the accelerator object for the x direction (may be NULL)
 * @param ya the accelerator object for the y direction (may be NULL)
 */
double interp2d_eval_deriv_xx(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the second y derivative of the interpolating function with the given data.
 * 
 * @param interp the interpolation object, previously initialized
 * @param xarr the x coordinates of the data, of length xsize
 * @param yarr the y coordinates of the data, of length ysize
 * @param zarr the z coordinates of the data, of length xsize*ysize
 * @param xsize the length of the array xa
 * @param ysize the length of the array ya
 * @param xa the accelerator object for the x direction (may be NULL)
 * @param ya the accelerator object for the y direction (may be NULL)
 */
double interp2d_eval_deriv_yy(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the cross derivative of the interpolating function with the given data.
 * 
 * @param interp the interpolation object, previously initialized
 * @param xarr the x coordinates of the data, of length xsize
 * @param yarr the y coordinates of the data, of length ysize
 * @param zarr the z coordinates of the data, of length xsize*ysize
 * @param xsize the length of the array xa
 * @param ysize the length of the array ya
 * @param xa the accelerator object for the x direction (may be NULL)
 * @param ya the accelerator object for the y direction (may be NULL)
 */
double interp2d_eval_deriv_xy(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Compute the index into a 1D array for a given pair of x,y indices.
 * This implements row-major ordering.
 */
#define INDEX_2D(xi, yi, xsize, ysize) ((yi) * (xsize) + (xi))

#ifdef __cplusplus
}
#endif

#endif
