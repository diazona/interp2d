`interp2d` generalizes the GSL interpolation routines to 2D interpolation. It tries to stick to the workflow and interface of the GSL 1D interpolation as closely as possible (even in cases where the GSL way may not quite make sense). The library includes implementations of bilinear and bicubic interpolation schemes.

The typical workflow is

1. create an interpolation object using `interp2d_alloc()`
2. initialize it using `interp2d_init()`
3. evaluate the interpolating function or its derivatives using `interp2d_eval()` or its counterparts, possibly many times
4. free the memory using `interp2d_free()`

If your code evaluates the interpolating function many times, it may benefit from using `gsl_interp_accel` objects, which are described [in the GSL documentation](http://www.gnu.org/software/gsl/manual/html_node/Index-Look_002dup-and-Acceleration.html#Index-Look_002dup-and-Acceleration).

`interp2d` uses CMake as its build system. Typically, to install a CMake-based project from source, from the directory containing `CMakeLists.txt`, run

    cmake .
    ./make
    ./make install

All files are licensed under the GPL version 3 or any later version at your option.
