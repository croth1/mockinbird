cmake_minimum_required(VERSION 2.8.12)
project(mockinbird CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(FindASan)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

add_executable(pcbam2mockinbird src/pcbam2mockinbird.cpp src/mockinbird_data.h)
add_executable(pcbam2regions src/pcbam2regions.cpp src/mockinbird_data.h)
