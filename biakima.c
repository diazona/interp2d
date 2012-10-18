/*
 * This file is part of interp2d, a GSL-compatible two-dimensional
 * interpolation library. <http://www.ellipsix.net/interp2d.html>
 *
 * Copyright 2012 Thomas Beutlich
 * Portions based on CALGO <http://www.netlib.org/toms/474>
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

#include <gsl/gsl_errno.h>
#include "interp2d.h"

#define dabs(x) (double)((x) >= 0 ? (x) : -(x))

static int itplbv_(size_t xsize, size_t ysize, const double xarr_[], const double yarr_[], const double zarr_[], double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double *z, double *z_x, double *z_y, double *z_xx, double *z_xy, double *z_yy) {
    double r__1;
    double za[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double zb[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double zab[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double p00[1] = {1};
    size_t ix, iy, jx, jy, jx1, jy1;
    double dx, dy;
    double zx[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double zy[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double zxy[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double w1[1] = {0}, w2[1] = {0};
    const size_t lxm1 = xsize - 1;
    const size_t lxm2 = lxm1 - 1;
    const size_t lxp1 = xsize + 1;
    const size_t lym1 = ysize - 1;
    const size_t lym2 = lym1 - 1;
    const size_t lyp1 = ysize + 1;
    double a3, b3;
    double x3, y3;
    double sw;
    double zx3b3;
    size_t jxml, jyml;
    double zx4b3, zy3a3, zy4a3;

#define a (zx + 1)
#define b (zx + 4)
#define c__ (zx + 7)
#define d__ (zy + 1)
#define e (zy + 13)
#define a1 (zx + 1)
#define b1 (zx + 1)
#define a2 (zx + 4)
#define a4 (zx + 7)
#define a5 (zx + 1)
#define b5 (zx + 1)
#define b2 (zy + 1)
#define b4 (zy + 13)
#define q0 (zx + 1)
#define q1 (zx + 4)
#define q2 (zx + 7)
#define q3 (zy + 1)
#define x2 (zx + 2)
#define y2 (zx + 13)
#define x4 (zx + 8)
#define x5 (zx + 11)
#define y4 (zy + 2)
#define y5 (zx + 14)
#define w4 (w2)
#define w3 (w1)
#define w5 (w1)
#define p01 (zy + 5)
#define p10 (zx + 5)
#define p11 (zxy + 5)
#define p02 (zx + 14)
#define p03 (zy + 4)
#define p12 (zy + 7)
#define p13 (zy + 8)
#define p20 (zy + 11)
#define p21 (zy + 14)
#define p22 (zxy + 1)
#define z33 (p00)
#define z23 (zy + 4)
#define z24 (zy + 7)
#define z32 (zy + 8)
#define z34 (zy + 11)
#define z35 (zy + 14)
#define z42 (zxy + 1)
#define z43 (zxy + 4)
#define p23 (zxy + 4)
#define z44 (zxy + 2)
#define p30 (zxy + 2)
#define z45 (zxy + 7)
#define p31 (zxy + 7)
#define z53 (zxy + 8)
#define p32 (zxy + 8)
#define z54 (zxy + 11)
#define p33 (zxy + 11)
#define wx2 (zxy + 13)
#define wy2 (w2)
#define wy3 (w1)
#define wx3 (zxy + 14)
#define z3a1 (za)
#define z3a2 (za + 1)
#define z3a3 (za + 2)
#define z3a4 (za + 3)
#define z3a5 (za + 4)
#define z4a1 (za + 5)
#define z4a2 (za + 6)
#define z4a3 (za + 7)
#define z4a4 (za + 8)
#define z4a5 (za + 9)
#define z3b1 (zb)
#define z3b2 (zb + 2)
#define z3b3 (zb + 4)
#define z3b4 (zb + 6)
#define z3b5 (zb + 8)
#define z4b1 (zb + 1)
#define z4b2 (zb + 3)
#define z4b3 (zb + 5)
#define z4b4 (zb + 7)
#define z4b5 (zb + 9)
#define zx33 (zx + 5)
#define zx43 (zx + 6)
#define zx34 (zx + 9)
#define zx44 (zx + 10)
#define zy33 (zy + 5)
#define zy43 (zy + 6)
#define zy34 (zy + 9)
#define zy44 (zy + 10)
#define a3sq (zx + 2)
#define b3sq (zy + 2)
#define jxm2 (jx1)
#define jym2 (jy1)
#define za2b2 (zab)
#define za3b2 (zab + 1)
#define za4b2 (zab + 2)
#define za2b3 (zab + 3)
#define za3b3 (zab + 4)
#define za4b3 (zab + 5)
#define za2b4 (zab + 6)
#define za3b4 (zab + 7)
#define za4b4 (zab + 8)
#define zxy33 (zxy + 5)
#define zxy43 (zxy + 6)
#define zxy34 (zxy + 9)
#define zxy44 (zxy + 10)

/* BIVARIATE INTERPOLATION */
/* THIS SUBROUTINE INTERPOLATES, FROM VALUES OF THE FUNCTION */
/* GIVEN AT INPUT GRID POINTS IN AN X-Y PLANE AND FOR A GIVEN */
/* SET OF POINTS IN THE PLANE, THE VALUES OF A SINGLE-VALUED */
/* BIVARIATE FUNCTION Z = Z(X,Y). */
/* THE METHOD IS BASED ON A PIECE-WISE FUNCTION COMPOSED OF */
/* A SET OF BICUBIC POLYNOMIALS IN X AND Y.  EACH POLYNOMIAL */
/* IS APPLICABLE TO A RECTANGLE OF THE INPUT GRID IN THE X-Y */
/* PLANE.  EACH POLYNOMIAL IS DETERMINED LOCALLY. */
/* THE INPUT PARAMETERS ARE */
/* IU  = LOGICAL UNIT NUMBER OF STANDARD OUTPUT UNIT */
/* LX  = NUMBER OF INPUT GRID POINTS IN THE X COORDINATE */
/*       (MUST BE 2 OR GREATER) */
/* LY  = NUMBER OF INPUT GRID POINTS IN THE Y COORDINATE */
/*       (MUST BE 2 OR GREATER) */
/* X   = ARRAY OF DIMENSION LX STORING THE X COORDINATES */
/*       OF INPUT GRID POINTS (IN ASCENDING ORDER) */
/* Y   = ARRAY OF DIMENSION LY STORING THE Y COORDINATES */
/*       OF INPUT GRID POINTS (IN ASCENDING ORDER) */
/* Z   = DOUBLY-DIMENSIONED ARRAY OF DIMENSION (LX,LY) */
/*       STORING THE VALUES OF THE FUNCTION (Z VALUES) */
/*       AT INPUT GRID POINTS */
/* N   = NUMBER OF POINTS AT WHICH INTERPOLATION OF THE */
/*       Z VALUE IS DESIRED (MUST BE 1 OR GREATER) */
/* U   = ARRAY OF DIMENSION N STORING THE X COORDINATES */
/*       OF DESIRED POINTS */
/* V   = ARRAY OF DIMENSION N STORING THE Y COORDINATES */
/*       OF DESIRED POINTS */
/* THE OUTPUT PARAMETER IS */
/* W   = ARRAY OF DIMENSION N WHERE THE INTERPOLATED Z */
/*       VALUES AT DESIRED POINTS ARE TO BE DISPLAYED */
/*       SOME VARIABLES INTERNALLY USED ARE */
/* ZA  = DIVIDED DIFFERENCE OF Z WITH RESPECT TO X */
/* ZB  = DIVIDED DIFFERENCE OF Z WITH RESPECT TO Y */
/* ZAB = SECOND ORDER DIVIDED DIFFERENCE OF Z WITH */
/*       RESPECT TO X AND Y */
/* ZX  = PARTIAL DERIVATIVE OF Z WITH RESPECT TO X */
/* ZY  = PARTIAL DERIVATIVE OF Z WITH RESPECT TO Y */
/* ZXY = SECOND ORDER PARTIAL DERIVATIVE OF Z WITH */
/*       RESPECT TO X AND Y */
/* DECLARATION STATEMENTS */
/* PRELIMINARY PROCESSING */
/* SETTING OF SOME INPUT PARAMETERS TO LOCAL VARIABLES */
    /* Parameter adjustments */
    const double *xarr = xarr_ - 1;
    const double *yarr = yarr_ - 1;
    const double *zarr = zarr_ - (1 + xsize);

/* ERROR CHECK */
    if (lxm2 < 0) {
        goto L710;
    }
    if (lym2 < 0) {
        goto L720;
    }

    // First compute the indices into the data arrays where we are interpolating
    if (xa != NULL) {
        ix = gsl_interp_accel_find(xa, xarr_, xsize, x);
    }
    else {
        ix = gsl_interp_bsearch(xarr_, x, 0, xsize - 1);
    }
    if (ya != NULL) {
        iy = gsl_interp_accel_find(ya, yarr_, ysize, y);
    }
    else {
        iy = gsl_interp_bsearch(yarr_, y, 0, ysize - 1);
    }
    ix += 2;
    iy += 2;

/* ROUTINES TO PICK UP NECESSARY X, Y, AND Z VALUES, TO */
/* COMPUTE THE ZA, ZB, AND ZAB VALUES, AND TO ESTIMATE THEM */
/* WHEN NECESSARY */
    jx = ix;
    if (jx == lxp1) {
        jx = xsize;
    }
    jy = iy;
    if (jy == lyp1) {
        jy = ysize;
    }
    jxm2 = jx - 2;
    jxml = jx - xsize;
    jym2 = jy - 2;
    jyml = jy - ysize;
/* IN THE CORE AREA, I.E., IN THE RECTANGLE THAT CONTAINS */
/* THE DESIRED POINT */
    x3 = xarr[jx - 1];
    *x4 = xarr[jx];
    a3 = 1. / (*x4 - x3);
    y3 = yarr[jy - 1];
    *y4 = yarr[jy];
    b3 = 1. / (*y4 - y3);
    *z33 = zarr[INDEX_2D(jx - 1, jy - 1, xsize, ysize)];
    *z43 = zarr[INDEX_2D(jx, jy - 1, xsize, ysize)];
    *z34 = zarr[INDEX_2D(jx - 1, jy, xsize, ysize)];
    *z44 = zarr[INDEX_2D(jx, jy, xsize, ysize)];
    *z3a3 = (*z43 - *z33) * a3;
    *z4a3 = (*z44 - *z34) * a3;
    *z3b3 = (*z34 - *z33) * b3;
    *z4b3 = (*z44 - *z43) * b3;
    *za3b3 = (*z4b3 - *z3b3) * a3;
/* IN THE X DIRECTION */
    if (lxm2 == 0) {
        goto L230;
    }
    if (jxm2 == 0) {
        goto L170;
    }
    *x2 = xarr[jx - 2];
    *a2 = 1. / (x3 - *x2);
    *z23 = zarr[INDEX_2D(jx - 2, jy - 1, xsize, ysize)];
    *z24 = zarr[INDEX_2D(jx - 2, jy, xsize, ysize)];
    *z3a2 = (*z33 - *z23) * *a2;
    *z4a2 = (*z34 - *z24) * *a2;
    if (jxml == 0) {
        goto L180;
    }
L170:
    *x5 = xarr[jx + 1];
    *a4 = 1. / (*x5 - *x4);
    *z53 = zarr[INDEX_2D(jx + 1, jy - 1, xsize, ysize)];
    *z54 = zarr[INDEX_2D(jx + 1, jy, xsize, ysize)];
    *z3a4 = (*z53 - *z43) * *a4;
    *z4a4 = (*z54 - *z44) * *a4;
    if (jxm2 != 0) {
        goto L190;
    }
    *z3a2 = *z3a3 + *z3a3 - *z3a4;
    *z4a2 = *z4a3 + *z4a3 - *z4a4;
    goto L190;
L180:
    *z3a4 = *z3a3 + *z3a3 - *z3a2;
    *z4a4 = *z4a3 + *z4a3 - *z4a2;
L190:
    *za2b3 = (*z4a2 - *z3a2) * b3;
    *za4b3 = (*z4a4 - *z3a4) * b3;
    if (jx <= 3) {
        goto L200;
    }
    *a1 = 1. / (*x2 - xarr[jx - 3]);
    *z3a1 = (*z23 - zarr[INDEX_2D(jx - 3, jy - 1, xsize, ysize)]) * *a1;
    *z4a1 = (*z24 - zarr[INDEX_2D(jx - 3 , jy, xsize, ysize)]) * *a1;
    goto L210;
L200:
    *z3a1 = *z3a2 + *z3a2 - *z3a3;
    *z4a1 = *z4a2 + *z4a2 - *z4a3;
L210:
    if (jx >= lxm1) {
        goto L220;
    }
    *a5 = 1. / (xarr[jx + 2] - *x5);
    *z3a5 = (zarr[INDEX_2D(jx + 2, jy - 1, xsize, ysize)] - *z53) * *a5;
    *z4a5 = (zarr[INDEX_2D(jx + 2, jy, xsize, ysize)] - *z54) * *a5;
    goto L240;
L220:
    *z3a5 = *z3a4 + *z3a4 - *z3a3;
    *z4a5 = *z4a4 + *z4a4 - *z4a3;
    goto L240;
L230:
    *z3a2 = *z3a3;
    *z4a2 = *z4a3;
    goto L180;
/* IN THE Y DIRECTION */
L240:
    if (lym2 == 0) {
        goto L310;
    }
    if (jym2 == 0) {
        goto L250;
    }
    *y2 = yarr[jy - 2];
    *b2 = 1. / (y3 - *y2);
    *z32 = zarr[INDEX_2D(jx - 1, jy - 2, xsize, ysize)];
    *z42 = zarr[INDEX_2D(jx, jy - 2, xsize, ysize)];
    *z3b2 = (*z33 - *z32) * *b2;
    *z4b2 = (*z43 - *z42) * *b2;
    if (jyml == 0) {
        goto L260;
    }
L250:
    *y5 = yarr[jy + 1];
    *b4 = 1. / (*y5 - *y4);
    *z35 = zarr[INDEX_2D(jx - 1, jy + 1, xsize, ysize)];
    *z45 = zarr[INDEX_2D(jx, jy + 1, xsize, ysize)];
    *z3b4 = (*z35 - *z34) * *b4;
    *z4b4 = (*z45 - *z44) * *b4;
    if (jym2 != 0) {
        goto L270;
    }
    *z3b2 = *z3b3 + *z3b3 - *z3b4;
    *z4b2 = *z4b3 + *z4b3 - *z4b4;
    goto L270;
L260:
    *z3b4 = *z3b3 + *z3b3 - *z3b2;
    *z4b4 = *z4b3 + *z4b3 - *z4b2;
L270:
    *za3b2 = (*z4b2 - *z3b2) * a3;
    *za3b4 = (*z4b4 - *z3b4) * a3;
    if (jy <= 3) {
        goto L280;
    }
    *b1 = 1. / (*y2 - yarr[jy - 3]);
    *z3b1 = (*z32 - zarr[INDEX_2D(jx - 1, jy - 3, xsize, ysize)]) * *b1;
    *z4b1 = (*z42 - zarr[INDEX_2D(jx, jy - 3, xsize, ysize)]) * *b1;
    goto L290;
L280:
    *z3b1 = *z3b2 + *z3b2 - *z3b3;
    *z4b1 = *z4b2 + *z4b2 - *z4b3;
L290:
    if (jy >= lym1) {
        goto L300;
    }
    *b5 = 1. / (yarr[jy + 2] - *y5);
    *z3b5 = (zarr[INDEX_2D(jx - 1, jy + 2, xsize, ysize)] - *z35) * *b5;
    *z4b5 = (zarr[INDEX_2D(jx, jy + 2, xsize, ysize)] - *z45) * *b5;
    goto L320;
L300:
    *z3b5 = *z3b4 + *z3b4 - *z3b3;
    *z4b5 = *z4b4 + *z4b4 - *z4b3;
    goto L320;
L310:
    *z3b2 = *z3b3;
    *z4b2 = *z4b3;
    goto L260;
/* IN THE DIAGONAL DIRECTIONS */
L320:
    if (lxm2 == 0) {
        goto L400;
    }
    if (lym2 == 0) {
        goto L410;
    }
    if (jxml == 0) {
        goto L350;
    }
    if (jym2 == 0) {
        goto L330;
    }
    *za4b2 = ((*z53 - zarr[INDEX_2D(jx + 1, jy - 2, xsize, ysize)]) * *b2 - *z4b2) * *a4;
    if (jyml == 0) {
        goto L340;
    }
L330:
    *za4b4 = ((zarr[INDEX_2D(jx + 1, jy + 1, xsize, ysize)] - *z54) * *b4 - *z4b4) * *a4;
    if (jym2 != 0) {
        goto L380;
    }
    *za4b2 = *za4b3 + *za4b3 - *za4b4;
    goto L380;
L340:
    *za4b4 = *za4b3 + *za4b3 - *za4b2;
    goto L380;
L350:
    if (jym2 == 0) {
        goto L360;
    }
    *za2b2 = (*z3b2 - (*z23 - zarr[INDEX_2D(jx - 2, jy - 2, xsize, ysize)]) * *b2) * *a2;
    if (jyml == 0) {
        goto L370;
    }
L360:
    *za2b4 = (*z3b4 - (zarr[INDEX_2D(jx - 2, jy + 1, xsize, ysize)] - *z24) * *b4) * *a2;
    if (jym2 != 0) {
        goto L390;
    }
    *za2b2 = *za2b3 + *za2b3 - *za2b4;
    goto L390;
L370:
    *za2b4 = *za2b3 + *za2b3 - *za2b2;
    goto L390;
L380:
    if (jxm2 != 0) {
        goto L350;
    }
    *za2b2 = *za3b2 + *za3b2 - *za4b2;
    *za2b4 = *za3b4 + *za3b4 - *za4b4;
    goto L420;
L390:
    if (jxml != 0) {
        goto L420;
    }
    *za4b2 = *za3b2 + *za3b2 - *za2b2;
    *za4b4 = *za3b4 + *za3b4 - *za2b4;
    goto L420;
L400:
    *za2b2 = *za3b2;
    *za4b2 = *za3b2;
    *za2b4 = *za3b4;
    *za4b4 = *za3b4;
    goto L420;
L410:
    *za2b2 = *za2b3;
    *za2b4 = *za2b3;
    *za4b2 = *za4b3;
    *za4b4 = *za4b3;
/* NUMERICAL DIFFERENTIATION   ---   TO DETERMINE PARTIAL */
/* DERIVATIVES ZX, ZY, AND ZXY AS WEIGHTED MEANS OF DIVIDED */
/* DIFFERENCES ZA, ZB, AND ZAB, RESPECTIVELY */
L420:
    for (jy = 2; jy <= 3; ++(jy)) {
        for (jx = 2; jx <= 3; ++(jx)) {
            *w2 = (r__1 = za[jx + 2 + (jy - 1) * 5 - 6] - za[jx + 1 + (jy - 1) * 5 - 6], dabs(r__1));
            *w3 = (r__1 = za[jx + (jy - 1) * 5 - 6] - za[jx - 1 + (jy - 1) * 5 - 6], dabs(r__1));
            sw = *w2 + *w3;
            if (sw == 0.) {
                goto L430;
            }
            *wx2 = *w2 / sw;
            *wx3 = *w3 / sw;
            goto L440;
L430:
            *wx2 = .5f;
            *wx3 = .5f;
L440:
            zx[jx + (jy << 2) - 5] = *wx2 * za[jx + (jy - 1) * 5 - 6] + *wx3 * za[jx + 1 + (jy - 1) * 5 - 6];
            *w2 = (r__1 = zb[jx - 1 + ((jy + 2) << 1) - 3] - zb[jx - 1 + ((jy + 1) << 1) - 3], dabs(r__1));
            *w3 = (r__1 = zb[jx - 1 + (jy << 1) - 3] - zb[jx - 1 + ((jy - 1) << 1) - 3], dabs(r__1));
            sw = *w2 + *w3;
            if (sw == 0.) {
                goto L450;
            }
            *wy2 = *w2 / sw;
            *wy3 = *w3 / sw;
            goto L460;
L450:
            *wy2 = .5f;
            *wy3 = .5f;
L460:
            zy[jx + (jy << 2) - 5] = *wy2 * zb[jx - 1 + (jy << 1) - 3] + *wy3 * zb[jx - 1 + ((jy + 1) << 1) - 3];
            zxy[jx + (jy << 2) - 5] = *wy2 * (*wx2 * zab[jx - 1 + (jy - 1) * 3 - 4] + *wx3 * zab[jx + (jy - 1) * 3 - 4]) + *wy3 * (*wx2 * zab[jx - 1 + jy * 3 - 4] + *wx3 * zab[jx + jy * 3 - 4]);
        }
    }
    if (ix == lxp1) {
        goto L530;
    }
    if (ix != 1) {
        goto L590;
    }
    *w2 = *a4 * (a3 * 3. + *a4);
    *w1 = a3 * 2. * (a3 - *a4) + *w2;
    for (jy = 2; jy <= 3; ++(jy)) {
        zx[(jy << 2) - 4] = (*w1 * za[(jy - 1) * 5 - 5] + *w2 * za[(jy - 1) * 5 - 4]) / (*w1 + *w2);
        zy[(jy << 2) - 4] = zy[(jy << 2) - 3] + zy[(jy << 2) - 3] - zy[(jy << 2) - 2];
        zxy[(jy << 2) - 4] = zxy[(jy << 2) - 3] + zxy[(jy << 2) - 3] - zxy[(jy << 2) - 2];
        for (jx1 = 2; jx1 <= 3; ++(jx1)) {
            jx = 5 - jx1;
            zx[jx + (jy << 2) - 5] = zx[jx - 1 + (jy << 2) - 5];
            zy[jx + (jy << 2) - 5] = zy[jx - 1 + (jy << 2) - 5];
            zxy[jx + (jy << 2) - 5] = zxy[jx - 1 + (jy << 2) - 5];
        }
    }
    x3 -= 1. / *a4;
    *z33 -= *z3a2 / *a4;
    for (jy = 1; jy <= 5; ++(jy)) {
        zb[(jy << 1) - 1] = zb[(jy << 1) - 2];
    }
    for (jy = 2; jy <= 4; ++(jy)) {
        zb[(jy << 1) - 2] -= zab[(jy - 1) * 3 - 3] / *a4;
    }
    a3 = *a4;
    jx = 1;
    goto L570;
L530:
    *w4 = *a2 * (a3 * 3. + *a2);
    *w5 = a3 * 2. * (a3 - *a2) + *w4;
    for (jy = 2; jy <= 3; ++(jy)) {
        zx[(jy << 2) - 1] = (*w4 * za[(jy - 1) * 5 - 2] + *w5 * za[(jy - 1) * 5 - 1]) / (*w4 + *w5);
        zy[(jy << 2) - 1] = zy[(jy << 2) - 2] + zy[(jy << 2) - 2] - zy[(jy << 2) - 3];
        zxy[(jy << 2) - 1] = zxy[(jy << 2) - 2] + zxy[(jy << 2) - 2] - zxy[(jy << 2) - 3];
        for (jx = 2; jx <= 3; ++(jx)) {
            zx[jx + (jy << 2) - 5] = zx[jx + 1 + (jy << 2) - 5];
            zy[jx + (jy << 2) - 5] = zy[jx + 1 + (jy << 2) - 5];
            zxy[jx + (jy << 2) - 5] = zxy[jx + 1 + (jy << 2) - 5];
        }
    }
    x3 = *x4;
    *z33 = *z43;
    for (jy = 1; jy <= 5; ++(jy)) {
        zb[(jy << 1) - 2] = zb[(jy << 1) - 1];
    }
    a3 = *a2;
    jx = 3;
L570:
    za[2] = za[jx];
    for (jy = 1; jy <= 3; ++(jy)) {
        zab[jy * 3 - 2] = zab[jx + jy * 3 - 4];
    }
/* WHEN (V(K).LT.Y(1)).OR.(V(K).GT.Y(LY)) */
L590:
    if (iy == lyp1) {
        goto L630;
    }
    if (iy != 1) {
        goto L680;
    }
    *w2 = *b4 * (b3 * 3. + *b4);
    *w1 = b3 * 2. * (b3 - *b4) + *w2;
    for (jx = 2; jx <= 3; ++(jx)) {
        if (jx == 3 && ix == lxp1) {
        goto L600;
        }
        if (jx == 2 && ix == 1) {
        goto L600;
        }
        zy[jx - 1] = (*w1 * zb[jx - 2] + *w2 * zb[jx]) / (*w1 + *w2);
        zx[jx - 1] = zx[jx + 3] + zx[jx + 3] - zx[jx + 7];
        zxy[jx - 1] = zxy[jx + 3] + zxy[jx + 3] - zxy[jx + 7];
L600:
        for (jy1 = 2; jy1 <= 3; ++(jy1)) {
        jy = 5 - jy1;
        zy[jx + (jy << 2) - 5] = zy[jx + ((jy - 1) << 2) - 5];
        zx[jx + (jy << 2) - 5] = zx[jx + ((jy - 1) << 2) - 5];
        zxy[jx + (jy << 2) - 5] = zxy[jx + ((jy - 1) << 2) - 5];
        }
    }
    y3 -= 1. / *b4;
    *z33 -= *z3b2 / *b4;
    *z3a3 -= *za3b2 / *b4;
    *z3b3 = *z3b2;
    *za3b3 = *za3b2;
    b3 = *b4;
    goto L670;
L630:
    *w4 = *b2 * (b3 * 3. + *b2);
    *w5 = b3 * 2. * (b3 - *b2) + *w4;
    for (jx = 2; jx <= 3; ++(jx)) {
        if (jx == 3 && ix == lxp1) {
            goto L640;
        }
        if (jx == 2 && ix == 1) {
            goto L640;
        }
        zy[jx + 11] = (*w4 * zb[jx + 4] + *w5 * zb[jx + 6]) / (*w4 + *w5);
        zx[jx + 11] = zx[jx + 7] + zx[jx + 7] - zx[jx + 3];
        zxy[jx + 11] = zxy[jx + 7] + zxy[jx + 7] - zxy[jx + 3];
L640:
        for (jy = 2; jy <= 3; ++(jy)) {
            zy[jx + (jy << 2) - 5] = zy[jx + ((jy + 1) << 2) - 5];
            zx[jx + (jy << 2) - 5] = zx[jx + ((jy + 1) << 2) - 5];
            zxy[jx + (jy << 2) - 5] = zxy[jx + ((jy + 1) << 2) - 5];
        }
    }
    y3 = *y4;
    *z33 += *z3b3 / b3;
    *z3a3 += *za3b3 / b3;
    *z3b3 = *z3b4;
    *za3b3 = *za3b4;
    b3 = *b2;
L670:
    if (ix != 1 && ix != lxp1) {
        goto L680;
    }
    jx = ix / lxp1 + 2;
    jx1 = 5 - jx;
    jy = iy / lyp1 + 2;
    jy1 = 5 - jy;
    zx[jx + (jy << 2) - 5] = zx[jx1 + (jy << 2) - 5] + zx[jx + (jy1 << 2) - 5] - zx[jx1 + (jy1 << 2) - 5];
    zy[jx + (jy << 2) - 5] = zy[jx1 + (jy << 2) - 5] + zy[jx + (jy1 << 2) - 5] - zy[jx1 + (jy1 << 2) - 5];
    zxy[jx + (jy << 2) - 5] = zxy[jx1 + (jy << 2) - 5] + zxy[jx + (jy1 << 2) - 5] - zxy[jx1 + (jy1 << 2) - 5];
/* DETERMINATION OF THE COEFFICIENTS OF THE POLYNOMIAL */
L680:
    zx3b3 = (*zx34 - *zx33) * b3;
    zx4b3 = (*zx44 - *zx43) * b3;
    zy3a3 = (*zy43 - *zy33) * a3;
    zy4a3 = (*zy44 - *zy34) * a3;
    *a = *za3b3 - zx3b3 - zy3a3 + *zxy33;
    *b = zx4b3 - zx3b3 - *zxy43 + *zxy33;
    *c__ = zy4a3 - zy3a3 - *zxy34 + *zxy33;
    *d__ = *zxy44 - *zxy43 - *zxy34 + *zxy33;
    *e = *a + *a - *b - *c__;
    *a3sq = a3 * a3;
    *b3sq = b3 * b3;
    *p02 = ((*z3b3 - *zy33) * 2. + *z3b3 - *zy34) * b3;
    *p03 = (*z3b3 * -2. + *zy34 + *zy33) * *b3sq;
    *p12 = ((zx3b3 - *zxy33) * 2. + zx3b3 - *zxy34) * b3;
    *p13 = (zx3b3 * -2. + *zxy34 + *zxy33) * *b3sq;
    *p20 = ((*z3a3 - *zx33) * 2. + *z3a3 - *zx43) * a3;
    *p21 = ((zy3a3 - *zxy33) * 2. + zy3a3 - *zxy43) * a3;
    *p22 = ((*a + *e) * 3. + *d__) * a3 * b3;
    *p23 = (*e * -3. - *b - *d__) * a3 * *b3sq;
    *p30 = (*z3a3 * -2. + *zx43 + *zx33) * *a3sq;
    *p31 = (zy3a3 * -2. + *zxy43 + *zxy33) * *a3sq;
    *p32 = (*e * -3. - *c__ - *d__) * b3 * *a3sq;
    *p33 = (*d__ + *e + *e) * *a3sq * *b3sq;
/* COMPUTATION OF THE POLYNOMIAL */
    dy = y - y3;
    *q0 = *p00 + dy * (*p01 + dy * (*p02 + dy * *p03));
    *q1 = *p10 + dy * (*p11 + dy * (*p12 + dy * *p13));
    *q2 = *p20 + dy * (*p21 + dy * (*p22 + dy * *p23));
    *q3 = *p30 + dy * (*p31 + dy * (*p32 + dy * *p33));
    dx = x - x3;
    if (z) {
        *z = *q0 + dx * (*q1 + dx * (*q2 + dx * *q3));
    }
    if (z_x) {
        *z_x = *q1 + dx * (2. * *q2 + 3. * dx * *q3);
    }
    if (z_y || z_xy) {
        double q1d = *p11 + dy * (2. * *p12 + 3. * dy * *p13);
        double q2d = *p21 + dy * (2. * *p22 + 3. * dy * *p23);
        double q3d = *p31 + dy * (2. * *p32 + 3. * dy * *p33);
        if (zy) {
            double q0d = *p01 + dy * (2. * *p02 + 3. * dy * *p03);
            *z_y = q0d + dx * (q1d + dx * (q2d + dx * q3d));
        }
        if (z_xy) {
            *z_xy = q1d + dx * (2. * q2d + 3. * dx * q3d);
        }
    }
    if (z_xx) {
        *z_xx = 2 * (*q2 + 3. * dx * *q3);
    }
    if (z_yy) {
        double q0dd = 2. * *p02 + 6. * dy * *p03;
        double q1dd = 2. * *p12 + 6. * dy * *p13;
        double q2dd = 2. * *p22 + 6. * dy * *p23;
        double q3dd = 2. * *p32 + 6. * dy * *p33;
        *z_yy = q0dd + dx * (q1dd + dx * (q2dd + dx * q3dd));
    }
/* NORMAL EXIT */
    return GSL_SUCCESS;
/* ERROR EXIT */
L710:
    GSL_ERROR_NULL("number of x values <= 1", GSL_EINVAL);
L720:
    GSL_ERROR_NULL("number of y values <= 1", GSL_EINVAL);
}

#undef zxy44
#undef zxy34
#undef zxy43
#undef zxy33
#undef za4b4
#undef za3b4
#undef za2b4
#undef za4b3
#undef za3b3
#undef za2b3
#undef za4b2
#undef za3b2
#undef za2b2
#undef jym2
#undef jxm2
#undef b3sq
#undef a3sq
#undef zy44
#undef zy34
#undef zy43
#undef zy33
#undef zx44
#undef zx34
#undef zx43
#undef zx33
#undef z4b5
#undef z4b4
#undef z4b3
#undef z4b2
#undef z4b1
#undef z3b5
#undef z3b4
#undef z3b3
#undef z3b2
#undef z3b1
#undef z4a5
#undef z4a4
#undef z4a3
#undef z4a2
#undef z4a1
#undef z3a5
#undef z3a4
#undef z3a3
#undef z3a2
#undef z3a1
#undef wx3
#undef wy3
#undef wy2
#undef wx2
#undef p33
#undef z54
#undef p32
#undef z53
#undef p31
#undef z45
#undef p30
#undef z44
#undef p23
#undef z43
#undef z42
#undef z35
#undef z34
#undef z32
#undef z24
#undef z23
#undef z33
#undef p22
#undef p21
#undef p20
#undef p13
#undef p12
#undef p03
#undef p02
#undef p11
#undef p10
#undef p01
#undef w5
#undef w3
#undef w4
#undef y5
#undef y4
#undef x5
#undef x4
#undef y2
#undef x2
#undef q3
#undef q2
#undef q1
#undef q0
#undef b4
#undef b2
#undef b5
#undef a5
#undef a4
#undef a2
#undef b1
#undef a1
#undef e
#undef d__
#undef c__
#undef b
#undef a

static int biakima_init(void* vstate, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize) {
    return GSL_SUCCESS;
}

static int biakima_eval(const void* vstate, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return itplbv_(xsize, ysize, xarr, yarr, zarr, x, y, xa, ya, z, NULL, NULL, NULL, NULL, NULL);
}

static int biakima_deriv_x(const void* vstate, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z_p) {
    return itplbv_(xsize, ysize, xarr, yarr, zarr, x, y, xa, ya, NULL, z_p, NULL, NULL, NULL, NULL);
}

static int biakima_deriv_y(const void* vstate, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z_p) {
    return itplbv_(xsize, ysize, xarr, yarr, zarr, x, y, xa, ya, NULL, NULL, z_p, NULL, NULL, NULL);
}

static int biakima_deriv_xx(const void* vstate, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z_pp) {
    return itplbv_(xsize, ysize, xarr, yarr, zarr, x, y, xa, ya, NULL, NULL, NULL, z_pp, NULL, NULL);
}

static int biakima_deriv_xy(const void* vstate, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z_pp) {
    return itplbv_(xsize, ysize, xarr, yarr, zarr, x, y, xa, ya, NULL, NULL, NULL, NULL, z_pp, NULL);
}

static int biakima_deriv_yy(const void* vstate, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z_pp) {
    return itplbv_(xsize, ysize, xarr, yarr, zarr, x, y, xa, ya, NULL, NULL, NULL, NULL, NULL, z_pp);
}

static const interp2d_type biakima_type = {
    "biakima",
    2,
    NULL,
    biakima_init,
    biakima_eval,
    biakima_deriv_x,
    biakima_deriv_y,
    biakima_deriv_xx,
    biakima_deriv_xy,
    biakima_deriv_yy,
    NULL
};

const interp2d_type* interp2d_biakima = &biakima_type;
