cmake_minimum_required(VERSION 3.0)

set(CMAKE_BUILD_TYPE Debug)

message("Cube Analysis Example -- ${VERSION}")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

add_executable(cubeExample.exe cubeExample.cxx)
target_link_libraries(cubeExample.exe LINK_PUBLIC CubeRecon::cuberecon_io)
install(TARGETS cubeExample.exe RUNTIME DESTINATION bin)

