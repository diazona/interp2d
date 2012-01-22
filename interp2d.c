#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "interp2d.h"

interp2d* interp2d_alloc(const interp2d_type* T, size_t size) {
    interp2d* interp;
    if (size < T->min_size) {
        GSL_ERROR_NULL("insufficient number of points for interpolation type", GSL_EINVAL);
    }
    interp = (interp2d*)malloc(sizeof(interp2d));
    if (interp == NULL) {
        GSL_ERROR_NULL("failed to allocate space for interp2d struct", GSL_ENOMEM);
    }
    interp->type = T;
    interp->size = size;
    if (interp->type->alloc == NULL) {
        interp->state = NULL;
        return interp;
    }
    interp->state = interp->type->alloc(size);
    if (interp->state == NULL) {
        free(interp);
        GSL_ERROR_NULL("failed to allocate space for interp2d state", GSL_ENOMEM);
    }
    return interp;
}

int interp2d_init(interp2d* interp, const double xarr[], const double yarr[], const double* zarr[], size_t size) {
    size_t i;
    if (size != interp->size) {
        GSL_ERROR("data must match size of interpolation object", GSL_EINVAL);
    }
    for (i = 1; i < size; i++) {
        if (xarr[i-1] >= xarr[i]) {
            GSL_ERROR("x values must be strictly increasing", GSL_EINVAL);
        }
        if (yarr[i-1] >= yarr[i]) {
            GSL_ERROR("y values must be strictly increasing", GSL_EINVAL);
        }
    }
    interp->xmin = xarr[0];
    interp->xmax = xarr[size - 1];
    interp->ymin = yarr[0];
    interp->ymax = yarr[size - 1];
    {
        int status = interp->type->init(interp->state, xarr, yarr, zarr, size);
        return status;
    }
}

void interp2d_free(interp2d* interp) {
    RETURN_IF_NULL(interp);
    if (interp->type->free) {
        interp->type->free(interp->state);
    }
    free(interp);
}

double interp2d_eval(const interp2d* interp, const double xarr[], const double yarr[], const double* zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    double z;
    int status;
    if (x < interp->xmin || x > interp->xmax) {
        GSL_ERROR_VAL("interpolation error", GSL_EDOM, GSL_NAN);
    }
    status = interp->type->eval(interp->state, xarr, yarr, zarr, interp->size, x, y, xa, ya, &z);
    DISCARD_STATUS(status);
    return z;
}

size_t interp2d_type_min_size(const interp2d_type* T) {
    return T->min_size;
}

size_t interp2d_min_size(const interp2d* interp) {
    return interp->type->min_size;
}

const char* interp2d_name(const interp2d* interp) {
    return interp->type->name;
}
