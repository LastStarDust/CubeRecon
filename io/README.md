# CubeRecon i/o library

This is the source for the cube reconstruction io class library.  This
is intended to be independent of the actual reconstruction code, and
it should be able to build outside of the cubeRecon package (or
without compiling the cube reconstruction code using the top level
CMakeList.txt file).  Programs that want to read the curbe
reconstruction results should include this library.
