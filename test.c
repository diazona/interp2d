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

int test_bilinear() {
    int status;
    size_t size = 4, test_size = 6;
    double xarr[4] = {0.0, 1.0, 2.0, 3.0};
    double yarr[4] = {0.0, 1.0, 2.0, 3.0};
    double zarr[16] = {1.0, 1.1, 1.2, 1.3,
                        1.1, 1.2, 1.3, 1.4,
                        1.2, 1.3, 1.4, 1.5,
                        1.3, 1.4, 1.5, 1.6};
    double xval[6] = {0.0, 0.5, 1.0, 1.5, 2.5, 3.0};
    double yval[6] = {0.0, 0.5, 1.0, 1.5, 2.5, 3.0};
    double zval[6] = {1.0, 1.1, 1.2, 1.3, 1.5, 1.6};
    status = test_interp2d(xarr, yarr, zarr, size, size, xval, yval, zval, test_size, interp2d_bilinear);
    gsl_test(status, "bilinear interpolation");
    return status;
}

int main(int argc, char** argv) {
    gsl_ieee_env_setup();
    argc = 0;
    argv = 0;
    test_bilinear();
    exit(gsl_test_summary());
}
