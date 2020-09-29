# CubeRecon

This is a reconstruction for cube based scintillation detectors.

# Running

The reconstruction is run in three stages.  The first stage translates
the input into the right format.  The second stage, and third stage
apply the reconstruction.

Files are translated into the right format by running
```cubeERepTranslate```.  This takes two arguments.  The first is the
name of an erep-sim output file, and the second is the name of the
translated CubeRecon file.

The second stage of the reconstruction is optional and will be run
later if the output doesn't exist.  However, it is often convenient to
apply this stage to simplify later debugging.  The second stage builds
the voxels out of the SiPM hits.  This will process a few events a
minute, so the time is not negligible.  The program is run as

```
cubeMakeHits3D input.root output.root
```

The final stage runs the actual reconstruction.  It will also build
the voxels if they don't exist.  It is run using the ```cubeRecon```
program.

```
cubeRecon input.root output.root
```

The reconstruction has not been optimized for efficiency, so you
should expect a huge amount of printed output, and that it will will
take about a minute per event.  Small events are much faster, but
large events are slow.

# Using the output

The output of the reconstruction is saved in the "CubeEvents" tree.
You can attach to the tree using a snippet of code similar to

```
std::unique_ptr<TFile> input(new TFile("input.root");
TTree* tree = std::dynamic_case<TTree*>(input->Get("CubeEvents"));
Cube::Event* event = NULL;
tree->SetBranchAddress("Event",&event);
tree->GetEntry(0);
```

The geometry used for the reconstruction is saved in a TGeoManager
object and can be accessed using

```
input->Get("CubeReconGeometry")
```

The "cuberecon_io" class library is necessary to access the events.
The classes are document in the relevant include files (...the source,
use the source...), but the important classes are:

- Cube::Event

- Cube::AlgorithmResult

- Cube::ReconTrack

- Cube::ReconCluster

## Add an example of a working analysis...

## Event display

The output of the reconstruction can be viewed using the ```cube-disp```
event display.  It is run with a single input file.  The display is
mostly aimed at debugging the reconstruction, so the GUI is not
static, but it will display most interesting features of an event.

# Building

The CubeRecon package uses CMake and the build is fairly standard.
There is a setup script located in the top directory that needs to be
able to find a working version of ROOT.  The setup script will make
sure that they can be located (using `thisroot.sh`).

```bash
. setup.sh
```

The package is built using cmake.  CMake can be run by hand, but there
is a script in the build directory that can be run using the
`cube-build` alias.  The build has been tested ROOT 6.22, but it will
probably build with any recent version ROOT.

# Requirements 

The only explicitly external requirement is that ROOT must be
available (and found by cmake).

If you are compiling ROOT and GEANT4 by hand, following will generally
work (changed for your versions):

```bash
tar xvzf root_v6.22.00.source.tar.gz
mkdir root-6.22.00-install root-6.22.00-build
cd root-6.22.00-build
cmake -DCMAKE_INSTALL_PREFIX=../root-6.22.00-install ../root-6.22.00
make
make install
```

### Expert Compilation

This is just a generic cmake build system, so everything can be done
by hand.  You need to make sure root is "in the path".  You can make
sure that ROOT is correctly setup using

``` bash
source thisroot.sh
```

Assuming that you have the source in "the-cube-recon-directory", will
build in "the-build-directory", and want to install in
"the-install-directory", the commands are:

```bash
cd the-build-directory
cmake -DCMAKE_INSTALL_PREFIX=the-install-directory the-cube-recon-directory
make
make install
```


