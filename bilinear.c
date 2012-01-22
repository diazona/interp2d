#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include "interp2d.h"

static int bilinear_init(void* state, const double xa[], const double ya[], const double za[], size_t size) {
    return GSL_SUCCESS;
}

static int bilinear_eval(const void* state, const double xarr[], const double yarr[], const double zarr[], size_t size, double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    double xmin, xmax, ymin, ymax, zminmin, zminmax, zmaxmin, zmaxmax;
    double dx, dy;
    double t, u;
    size_t xi, yi;
    if (xa != NULL) {
        xi = gsl_interp_accel_find(xa, xarr, size, x);
    }
    else {
        xi = gsl_interp_bsearch(xarr, x, 0, size - 1);
    }
    if (ya != NULL) {
        yi = gsl_interp_accel_find(ya, yarr, size, y);
    }
    else {
        yi = gsl_interp_bsearch(yarr, y, 0, size - 1);
    }
    xmin = xarr[xi];
    xmax = xarr[xi + 1];
    ymin = yarr[yi];
    ymax = yarr[yi + 1];
    zminmin = zarr[INDEX_2D(xi, yi, size, size)];
    zminmax = zarr[INDEX_2D(xi, yi + 1, size, size)];
    zmaxmin = zarr[INDEX_2D(xi + 1, yi, size, size)];
    zmaxmax = zarr[INDEX_2D(xi + 1, yi + 1, size, size)];
    dx = xmax - xmin;
    dy = ymax - ymin;
    t = (x - xmin)/dx;
    u = (y - ymin)/dy;
    *z = (1.-t)*(1.-u)*zminmin + t*(1.-u)*zmaxmin + (1.-t)*u*zminmax + t*u*zmaxmax;
    return GSL_SUCCESS;
}

static const interp2d_type bilinear_type = {
    "bilinear",
    2,
    NULL,
    &bilinear_init,
    &bilinear_eval
};

const interp2d_type* interp2d_bilinear = &bilinear_type;