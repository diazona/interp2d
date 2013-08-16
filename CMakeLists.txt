cmake_minimum_required(VERSION 2.6)
project(interp2d)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Code snippet from http://majewsky.wordpress.com/2010/08/14/tip-of-the-day-cmake-and-doxygen/
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
                    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

add_library(interp2d bilinear.c bicubic.c interp2d.c interp2d_spline.c)
install(TARGETS interp2d DESTINATION lib)
install(FILES interp2d.h interp2d_spline.h DESTINATION include)

add_executable(interp2dtest test.c)
target_link_libraries(interp2dtest m gsl interp2d)

enable_testing()
add_test(all interp2dtest)
