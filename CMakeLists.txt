cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(CMAKE_CXX_FLAGS "-O3 -lm -fPIC -fomit-frame-pointer -Wall -Wextra")

project(rct)

find_package(HDF5 REQUIRED COMPONENTS C CXX)

add_library(rct STATIC src/rct.cc src/random.cc src/vec.cc)
set_property(TARGET rct PROPERTY CXX_STANDARD 17)

add_executable(build src/build.cc)
target_link_libraries(build rct ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
set_property(TARGET build PROPERTY CXX_STANDARD 17)
include_directories(${HDF5_INCLUDE_DIRS})
