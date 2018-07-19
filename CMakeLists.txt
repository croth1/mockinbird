cmake_minimum_required(VERSION 2.8.12)
project(mockinbird CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

add_executable(iclipbam2mockinbird src/iclipbam2mockinbird.cpp)
