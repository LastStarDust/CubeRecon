#!/bin/sh
#
# Setup the CubeRecon directory for development (or simply running the
# reconstruction).  This makes sure that the ROOT environment
# variables are set (using thisroot.sh), and the edep-sim i/o library
# can be found by cmake.  This defines a couple of conveniences
# commands:
#
#  cube-build == Source ./build/cube-build.sh which will conveniently
#           run cmake/make/make install from any place so that it's
#           really easy to recompile.
#
#  cube-setup == Source this file.  You probably never have to use
#           this one.
#
# This setup script is not needed.  You can also do it by hand.  It's
# a usual cmake build, but you need to make sure root is
# "in the path".
#
# source thisroot.sh
# cd the-build-directory
# cmake -DCMAKE_INSTALL_PREFIX=the-install-directory the-cube-reco-directory 
# make
# make install


# Try to setup root.  ROOT installs thisroot.sh in the bin directory
# to setup the environment.  The first "thisroot.sh" in the path will
# define the root that is used.
. thisroot.sh >& /dev/null
if [ $? != 0 ]; then
    echo ROOT not available.
fi

# Find the root of the building area.
___cube_root() {
    COUNT=50
    while true; do
	if [ -e ./build -a -d ./build -a -e ./build/cube-build.sh ]; then
	    echo ${PWD}
	    return
	fi
	COUNT=$(expr ${COUNT} - 1)
	if [ ${COUNT} -lt 1 ]; then
	    echo invalid-directory
	    return
	fi
	cd ..
    done
}

export CUBE_ROOT
CUBE_ROOT=$(___cube_root)
unset -f ___cube_root

if [ ${CUBE_ROOT} = "invalid-directory" ]; then
    echo The CubeRecon setup.sh must be sourced in cube-reco directory.
    return
fi

___cube_target () {
    target="cube"
    compiler=gcc
    target="${target}-${compiler}-$(cc -dumpversion)-$(cc -dumpmachine)"
    echo $target
}

export CUBE_TARGET
CUBE_TARGET=$(___cube_target)
unset -f ___cube_target

___path_prepend () {
    ___path_remove $1 $2
    eval export $1="$2:\$$1"
}
___path_remove ()  {
    export $1=$(eval echo -n \$$1 | \
	awk -v RS=: -v ORS=: '$0 != "'$2'"' | \
	sed 's/:$//'); 
}

___path_prepend PATH ${CUBE_ROOT}/${CUBE_TARGET}/bin
___path_prepend LD_LIBRARY_PATH ${CUBE_ROOT}/${CUBE_TARGET}/lib

unset -f ___path_prepend
unset -f ___path_remove


alias cube-setup=". ${CUBE_ROOT}/setup.sh"

alias cube-build="${CUBE_ROOT}/build/cube-build.sh"

echo Defined cube-setup to re-setup the cube-reco package.
echo Defined cube-build to build the the cube-reco package.
