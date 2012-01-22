project(interp2d)

add_executable(interp2dtest main.cpp)
add_library(interp2d bilinear.c interp2d.c)

