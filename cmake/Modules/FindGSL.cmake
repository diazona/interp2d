# Tell CMake how to find GSL
# Once done this will define
#  GSL_FOUND         - true
#  GSL_INCLUDE_DIRS  - include path for GSL headers
#  GSL_LIBRARIES     - link path for GSL libraries

include(LibFindMacros)

libfind_pkg_check_modules(GSL_PKGCONF GSL)

find_path(GSL_INCLUDE_DIR
  NAMES gsl/gsl_errno.h gsl/gsl_integration.h gsl/gsl_interp.h gsl/gsl_math.h
        gsl/gsl_monte.h gsl/gsl_monte_miser.h gsl/gsl_monte_vegas.h gsl/gsl_qrng.h
        gsl/gsl_rng.h gsl/gsl_roots.h gsl/gsl_sf.h gsl/gsl_sf_bessel.h gsl/gsl_sys.h
  PATHS ${GSL_PKGCONF_INCLUDE_DIRS}
)

find_library(GSL_LIBRARY
  NAMES gsl
  PATHS ${GSL_PKGCONF_LIBRARY_DIRS}
)

find_library(GSL_CBLAS_LIBRARY
  NAMES gslcblas
  PATHS ${GSL_PKGCONF_LIBRARY_DIRS}
)

find_library(CMATH_LIBRARY
  NAMES m
)

set(GSL_PROCESS_INCLUDES GSL_INCLUDE_DIR GSL_INCLUDE_DIRS)
set(GSL_PROCESS_LIBS GSL_LIBRARY GSL_CBLAS_LIBRARY CMATH_LIBRARY GSL_LIBRARIES)
libfind_process(GSL)
