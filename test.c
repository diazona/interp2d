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

#include <math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include "interp2d.h"

int test_interp2d(const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize, const double xval[], const double yval[], const double zval[], size_t test_size, const interp2d_type* T) {
    gsl_interp_accel *xa, *ya;
    int status = 0;
    size_t i;

    xa = gsl_interp_accel_alloc();
    ya = gsl_interp_accel_alloc();
    interp2d* interp = interp2d_alloc(T, xsize, ysize);
    unsigned int min_size = interp2d_type_min_size(T);
    
    gsl_test_int(min_size, T->min_size, "interp2d_type_min_size on %s", interp2d_name(interp));
    interp2d_init(interp, xarr, yarr, zarr, xsize, ysize);
    for (i = 0; i < test_size; i++) {
        double x = xval[i];
        double y = yval[i];
        double z;
        z = interp2d_eval(interp, xarr, yarr, zarr, x, y, xa, ya);
        gsl_test_abs(z, zval[i], 1e-10, "%s %d", interp2d_name(interp), i);
        if (fabs(z - zval[i]) > 1e-10) {
            status++;
        }
    }
    gsl_interp_accel_free(xa);
    gsl_interp_accel_free(ya);
    interp2d_free(interp);
    return status;
}

int test_bilinear_symmetric() {
    int status;
    double xarr[] = {0.0, 1.0, 2.0, 3.0};
    double yarr[] = {0.0, 1.0, 2.0, 3.0};
    double zarr[] = {1.0, 1.1, 1.2, 1.3,
                     1.1, 1.2, 1.3, 1.4,
                     1.2, 1.3, 1.4, 1.5,
                     1.3, 1.4, 1.5, 1.6};
    double xval[] = {0.0, 0.5, 1.0, 1.5, 2.5, 3.0};
    double yval[] = {0.0, 0.5, 1.0, 1.5, 2.5, 3.0};
    double zval[] = {1.0, 1.1, 1.2, 1.3, 1.5, 1.6};
    size_t xsize = sizeof(xarr) / sizeof(xarr[0]);
    size_t ysize = sizeof(yarr) / sizeof(yarr[0]);
    size_t test_size = sizeof(xval) / sizeof(xval[0]);
    status = test_interp2d(xarr, yarr, zarr, xsize, ysize, xval, yval, zval, test_size, interp2d_bilinear);
    gsl_test(status, "bilinear interpolation with symmetric values");
    return status;
}

int test_bilinear_asymmetric_z() {
    int status;
    double xarr[] = {0.0, 1.0, 2.0, 3.0};
    double yarr[] = {0.0, 1.0, 2.0, 3.0};
    double zarr[] = {1.0, 1.3, 1.5, 1.6,
                     1.1, 1.4, 1.6, 1.9,
                     1.2, 1.5, 1.7, 2.2,
                     1.4, 1.7, 1.9, 2.3};
    double xval[] = {0.0, 0.5, 1.0, 1.5,  2.5,   3.0,  1.3954,    1.6476,       0.824957,  2.41108,  2.98619,   1.36485};
    double yval[] = {0.0, 0.5, 1.0, 1.5,  2.5,   3.0,  0.265371,  2.13849,      1.62114,   1.22198,  0.724681,  0.0596087};
    // results computed using Mathematica 8.0.1
    double zval[] = {1.0, 1.2, 1.4, 1.55, 2.025, 2.3,  1.2191513, 1.7242442248, 1.5067237, 1.626612, 1.6146423, 1.15436761};
    size_t xsize = sizeof(xarr) / sizeof(xarr[0]);
    size_t ysize = sizeof(yarr) / sizeof(yarr[0]);
    size_t test_size = sizeof(xval) / sizeof(xval[0]);
    status = test_interp2d(xarr, yarr, zarr, xsize, ysize, xval, yval, zval, test_size, interp2d_bilinear);
    gsl_test(status, "bilinear interpolation with asymmetric z values");
    return status;
}

int main(int argc, char** argv) {
    gsl_ieee_env_setup();
    argc = 0;
    argv = 0;
    test_bilinear_symmetric();
    test_bilinear_asymmetric_z();
    exit(gsl_test_summary());
}
