#ifndef __INTERP_2D_H__
#define __INTERP_2D_H__

#include <gsl/gsl_interp.h>

typedef struct {
    const char* name;
    unsigned int min_size;
    void* (*alloc)(size_t size);
    int (*init)(void*, const double xa[], const double ya[], const double* za[], size_t size);
    int (*eval)(const void*, const double xa[], const double ya[], const double* za[], size_t size, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z);
//     int (*eval_deriv) (const void*, const double xa[], const double ya[], const double za[][], size_t size, double x, double y, interp2d_accel*, double* z_p);
//     int (*eval_deriv2) (const void*, const double xa[], const double ya[], const double za[][], size_t size, double x, double y, interp2d_accel*, double* z_pp);
//     int (*integ) (const void*, const double xa[], const double ya[], const double za[][], size_t size, interp2d_accel*, double xa, double xb, double ya, double yb, double* result);
    void (*free)(void*);
} interp2d_type;

typedef struct {
    const interp2d_type* type;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    size_t size;
    void* state;
} interp2d;

GSL_VAR const interp2d_type* interp2d_bilinear;

interp2d* interp2d_alloc(const interp2d_type* T, size_t size);
int interp2d_init(interp2d* interp, const double xa[], const double ya[], const double* za[], size_t size);
void interp2d_free(interp2d* interp);

double interp2d_eval(const interp2d* interp, const double xarr[], const double yarr[], const double* zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

#endif
