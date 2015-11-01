`interp2d` generalizes the GSL interpolation routines to 2D interpolation. It tries to stick to the workflow and interface of the GSL 1D interpolation as closely as possible (even in cases where the GSL way may not quite make sense). The library includes implementations of bilinear and bicubic interpolation schemes.

As of GSL 2.0, the code from `interp2d` (bilinear and bicubic interpolation) is part of the GSL itself. There is no need to use this library with GSL 2.x.

Just like the GSL interpolation functions, there are two interfaces to the code. You can use the low-level interface defined in interp2d.h, which does not store the data arrays that define the function being interpolated (so you have to store them yourself and pass them to every function call), or you can use the high-level interface in interp2d_spline.h, which does store the data arrays in the interp2d_spline object.

The typical workflow is

1. create an interpolation object using `interp2d_alloc()` (low-level) or `interp2d_spline_alloc()` (high-level)
2. initialize it using `interp2d_init()` or `interp2d_spline_alloc()`
3. evaluate the interpolating function or its derivatives using `interp2d_eval()`/`interp2d_spline_eval()` or its counterparts, possibly many times
4. free the memory using `interp2d_free()` or `interp2d_spline_free()`

If your code evaluates the interpolating function many times, it may benefit from using `gsl_interp_accel` objects, which are described [in the GSL documentation](http://www.gnu.org/software/gsl/manual/html_node/Index-Look_002dup-and-Acceleration.html#Index-Look_002dup-and-Acceleration).

`interp2d` uses CMake as its build system. Typically, to install a CMake-based project from source, from the directory containing `CMakeLists.txt`, run

    cmake .
    ./make
    ./make install

All files are licensed under the GPL version 3 or any later version at your option.

Development of `interp2d` is hosted on the Github project page <https://github.com/diazona/interp2d>, which is where you should report bugs or offer suggestions for improvement if you would like to do so.
