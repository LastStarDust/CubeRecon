# Example for reading the CUBE_RECON output

This is an example of a very simple histogramming program for cube
recon.  It needs to have CubeRecon in the CMAKE_PREFIX_PATH so it can
find the cuberecon_io library, but other than that it can be copied
and used any place.  You can start using this by first setting up
ROOT, then CubeRecon, and finally this package

```
source thisroot.sh
cd CubeRecon
source setup.sh
cd example
source setup.sh
```

This subdirectory should not be built as part of the main CubeRecon
build (i.e. don't add it to ```../CMakeLists.txt```)
