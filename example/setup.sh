#!/bin/sh
#
# Setup the example directory for development (or simply running the
# simulation).  This makes sure that the ROOT environment variables
# are set (using thisroot.sh), and then defines a couple of
# conveniences commands:
#
#  example-build == Source ./build/example-build.sh which will conveniently
#           run cmake/make/make install from any place so that it's
#           really easy to recompile.
#
#  example-setup == Source this file.  You probably never have to use
#           this one.
#
# This setup script is not needed.  You can also do it by hand.  It's
# a usual cmake build, but you need to make sure root and geant are
# "in the path".
#
# source thisroot.sh
# cd the-build-directory
# cmake -DCMAKE_INSTALL_PREFIX=the-install-directory the-example-directory 
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
___example_root() {
    COUNT=50
    while true; do
	if [ -e ./build -a -d ./build -a -e ./build/example-build.sh ]; then
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

export EXAMPLE_ROOT
EXAMPLE_ROOT=$(___example_root)
unset -f ___example_root

if [ ${EXAMPLE_ROOT} = "invalid-directory" ]; then
    echo EXAMPLE setup.sh must be sourced in example directory.
    return
fi

___example_target () {
    target="example"
    compiler=gcc
    target="${target}-${compiler}-$(cc -dumpversion)-$(cc -dumpmachine)"
    echo $target
}

export EXAMPLE_TARGET
EXAMPLE_TARGET=$(___example_target)
unset -f ___example_target

___path_prepend () {
    ___path_remove $1 $2
    eval export $1="$2:\$$1"
}
___path_remove ()  {
    export $1=$(eval echo -n \$$1 | \
	awk -v RS=: -v ORS=: '$0 != "'$2'"' | \
	sed 's/:$//'); 
}

___path_prepend PATH ${EXAMPLE_ROOT}/${EXAMPLE_TARGET}/bin
___path_prepend LD_LIBRARY_PATH ${EXAMPLE_ROOT}/${EXAMPLE_TARGET}/lib

unset -f ___path_prepend
unset -f ___path_remove


alias example-setup=". ${EXAMPLE_ROOT}/setup.sh"

alias example-build="${EXAMPLE_ROOT}/build/example-build.sh"

echo Defined example-setup to re-setup the example package.
echo Defined example-build to build the the example package.
