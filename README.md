# sos_flow_experiments

This repository contains several examples that use SOS and/or TAU for runtime monitoring of parallel applications and workflows.

SOS has an optional dependency on EVPath. Several of the workflows use ADIOS.  For specific information on installing ADIOS and its dependencies (i.e. EVPath), see [https://www.olcf.ornl.gov/center-projects/adios/](https://www.olcf.ornl.gov/center-projects/adios/).  In the instructions below, the "CHAOS" path is the set of dependencies for ADIOS (EVPath, etc.)

# Building ADIOS for Pooky event extraction

When working with ADIOS applications, you should use this version of ADIOS (for now): [https://github.com/khuck/ADIOS] (https://github.com/khuck/ADIOS).  It includes support for extracting the dimensions at each adios_write call.

NOTE:
You will also need to use this version of TAU (for now): [http://www.nic.uoregon.edu/~khuck/tau2-git-pooky.tar.gz] (http://www.nic.uoregon.edu/~khuck/tau2-git-pooky.tar.gz).  

# Building SOS_flow

### To configure and install SOS (quick):

```
git clone https://github.com/cdwdirect/sos_flow.git
cd sos_flow
source hosts/linux/setenv.sh
./scripts/configure.sh -c
cd build-linux
make && make install
```

### To configure and build SOS (if quick doesn't work - modify as necessary to match your filesystem and compilers):

```
git clone https://github.com/cdwdirect/sos_flow.git
cd sos_flow
export CHAOS=$HOME/src/chaos/linux-gcc
export PKG_CONFIG_PATH=${CHAOS}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${CHAOS}/lib:${LD_LIBRARY_PATH}
export PATH=${CHAOS}/bin:${ADIOS_ROOT}/bin:${PATH}
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DSOSD_CLOUD_SYNC_WITH_EVPATH=TRUE -DSOSD_CLOUD_SYNC_WITH_MPI=FALSE -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx /gpfs/home/khuck/src/sos_flow/scripts/..
make && make install
```

# Building TAU and PDT

### Building PDT

One of the optional dependencies of TAU is PDT:

```
wget http://tau.uoregon.edu/pdt.tgz
tar -xvzf pdt.tgz
cd pdtoolkit-3.24
./configure -prefix=/usr/local/pdtoolkit/3.24
# this step may require sudo
make && make install
```

### Building TAU

To configure & build TAU, use this patched version of TAU: [http://www.nic.uoregon.edu/~khuck/tau2-git-latest.tar.gz](http://www.nic.uoregon.edu/~khuck/tau2-git-latest.tar.gz).  The paths to SOS and ADIOS are examples, please modify for your filesystem.

```
wget http://www.nic.uoregon.edu/~khuck/tau2-git-pooky.tar.gz
tar -xvzf http://www.nic.uoregon.edu/~khuck/tau2-git-pooky.tar.gz
cd tau2-git-latest
./configure -adios=/home/khuck/src/chaos/adios/ADIOS-gcc -sos=/home/khuck/src/sos_flow/build-linux -pdt=/usr/local/packages/pdt/3.23 -mpi -pthread
```

# Building and Running an example

The example of interest is a coupled application that uses ADIOS to exchange data from one MPI application to another.  TAU will intercept the MPI and ADIOS events, and send the performance data to SOS, which will aggregate the data in a database.  After the example executes, we can post-process the database to extract the MPI collective periodicity.

### Building the example

The example of interest is in the "Junmin" directory of this repo.

```
cd Junmin
./build-linux.sh
```

The key change to building xmain and reader2 is that the TAU libraries need to be added *before* the ADIOS libraries in the link line.  TAU will perform the measurements, and send the TAU data to the SOS aggregator network.  See ImpactTv1betaAdios/Makefile.tau and readerFull/compile.tau for details:

```
ADIOS = $(shell tau_cc.sh -tau:showlibs) $(shell ${ADIOSDIR}/bin/adios_config -l -f)
```

### Running the example

```
mkdir testdir
# copy the executables
cp ImpactTv1betaAdios/xmain testdir
cp readerFull/reader2 testdir
# copy the run scripts/inputs
cp tiny/ImpactT.in testdir
cp tiny/small.r testdir
cd testdir
# modify small.r as necessary
sbatch small.r
```

### Postprocessing

To extract the MPI collective periodicity from the SOS database, run a script in the SOS_flow repo (modify the path as necessary).  The script takes two arguments, the first is the name of the SOS database, and the second argument is the name of the application of interest (since the SOS database contains performance data from both applications):

```
python $HOME/src/sos_flow/src/soapy/compute-window-periodicity.py ./sosd.00000.db xmain
```

The script will open the database, perform queries and some analysis, and should eventually make output like this:

```
periodicity:  1.90189814568
compute duration, first arrival:  0.560921251774 , last arrival:  0.760286927223
last window: 
1499297094.66 to 1499297095.22
next 3 windows: 
1499297096.49 to 1499297097.05
1499297098.31 to 1499297098.87
1499297100.13 to 1499297100.69
Closing connection to database.
Done.
```

The script will also open a chart window that will look something like this:

![./figures/mpi_periodicity.png](./figures/mpi_periodicity.png)

The script will also (eventuallty) generate a YAML file with the MPI collective events, to generate an ADIOS Skel/Diesel proxy application (TBD).
