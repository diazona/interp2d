cmake_minimum_required(VERSION 2.6)
project(interp2d)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

add_library(interp2d bilinear.c interp2d.c)

add_executable(interp2dtest test.c)
target_link_libraries(interp2dtest m gsl interp2d)
