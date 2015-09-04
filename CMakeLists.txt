cmake_minimum_required(VERSION 2.6)
project(interp2d)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(LIBS ${LIBS} m)

# Code snippet from http://majewsky.wordpress.com/2010/08/14/tip-of-the-day-cmake-and-doxygen/
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
                    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

option(interp2d_STATIC_LIBRARY "Build interp2d as a static library" ON)
option(interp2d_SHARED_LIBRARY "Build interp2d as a shared library" OFF)

if (interp2d_STATIC_LIBRARY)
  add_library(interp2d_static bilinear.c bicubic.c interp2d.c interp2d_spline.c)
  set_target_properties(interp2d_static PROPERTIES OUTPUT_NAME interp2d)
  install(TARGETS interp2d_static DESTINATION lib)
endif (interp2d_STATIC_LIBRARY)
if (interp2d_SHARED_LIBRARY)
  add_library(interp2d_shared bilinear.c bicubic.c interp2d.c interp2d_spline.c)
  if (NOT interp2d_STATIC_LIBRARY)
    set_target_properties(interp2d_shared PROPERTIES OUTPUT_NAME interp2d)
  endif (NOT interp2d_STATIC_LIBRARY)
  install(TARGETS interp2d_shared DESTINATION lib)
endif (interp2d_SHARED_LIBRARY)

install(FILES interp2d.h interp2d_spline.h DESTINATION include)

# FindGSL.cmake does not give full paths to libraries like it's supposed
# to. (Blame gsl-config, which does not obey the CMake convention of
# returning full paths to libraries.) So it's necessary to explicitly
# add directories in which to search for the GSL libraries.
link_directories(${GSL_LIBRARY_DIRS})

add_executable(interp2dtest test.c)
if (interp2d_STATIC_LIBRARY)
  target_link_libraries(interp2dtest interp2d_static ${LIBS})
else (interp2d_STATIC_LIBRARY)
  target_link_libraries(interp2dtest interp2d_shared ${LIBS})
endif (interp2d_STATIC_LIBRARY)

enable_testing()
add_test(all interp2dtest)
