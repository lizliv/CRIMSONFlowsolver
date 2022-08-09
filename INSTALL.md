The CRIMSON flowsolver can be obtained from here https://github.com/abhilashreddym/CRIMSONFlowsolver.git. (main branch)

# Building CRIMSON flowsolver on Ubuntu
These are written for Ubuntu, but may be adapted for typical linux on cluster computers. Some of the dependencies that require sudo will be available through environment modules or similar. The rest would have to be built from source.

## Install dependencies from repo
```bash
sudo apt update
sudo apt install build-essential cmake git openmpi-bin libopenmpi-dev libboost-thread-dev libboost-filesystem-dev libboost-dev libreadline-dev libncurses-dev
```

The rest of the dependencies will have to be either built from source or obtained. In these instructions,  `~/temp` is used as the build folder
```bash
mkdir ~/temp
```
## Openblas

Download the source for the latest release from https://github.com/xianyi/OpenBLAS/releases. Below I have used the latest one as I write this.

```bash
cd ~/temp
mkdir openblas && cd $_ #make directory and change in to it
wget https://github.com/xianyi/OpenBLAS/releases/download/v0.3.20/OpenBLAS-0.3.20.tar.gz # download the source
tar xf OpenBLAS-0.3.20.tar.gz # unzip the archive
cd OpenBLAS-0.3.20
make -j4  # build openblas
make PREFIX=$HOME/.local/openblas install # install the library and include files in $HOME:/.local
ls ~/.local/openblas/*  # check that the files got installed properly
```

## PetSc
Now we will get and build PETSc.
```bash
cd ~/temp
git clone -b release https://gitlab.com/petsc/petsc.git petsc
```
Build PetSc
```bash
cd petsc
./configure  --prefix=$HOME/.local/petsc/ --with-memalign=64 COPTFLAGS='-O3 -march=native' FOPTFLAGS='-O3 -march=native'  CXXOPTFLAGS='-O3 -march=native' --with-x=0 --with-blaslapack-lib=$HOME/.local/openblas/lib/libopenblas.a --with-debugging=0
```

Once configure is complete, it will print a message with the command to build similar to below. run this command to begin the build. You can add an option `-j2` or `-j4` to make it go a bit faster but not more than that so as not to hog the login node.
**Dont copy the below command. Use the command that configure prints at the end**
```bash
make -j4 PETSC_DIR=$HOME/temp/petsc/petsc PETSC_ARCH=arch-linux-c-opt all
```

Once the build is complete you will see the command to install the files. **again dont copy the below command. use the one the make prints at the end**
```bash
make PETSC_DIR=$HOME/temp/petsc/petsc PETSC_ARCH=arch-linux-c-opt install
```

This will install the files at the prefix location `$HOME/.local/petsc`.
`make` will suggest that you do a check to verify everything is ok, which is recommended.

## VTK
Get the CRIMSON modified vtk from https://github.com/carthurs/VTK-5.8.0
```bash
cd ~/temp
git clone --depth=1 https://github.com/carthurs/VTK-5.8.0.git
cd VTK-5.8.0
cmake -S . -B . -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DVTK_USE_RENDERING=OFF -DCMAKE_INSTALL_PREFIX=$HOME/.local/vtk-crimson -DCMAKE_CXX_FLAGS=-Wno-error=narrowing
make -j4 && make install
```

## Python 2.7
CRIMSONFlowsolver needs libpython2.7 and associated `Python.h` header file. You could install it using the package manager, conda or build from source. . 

### Using Conda
Liz suggested the this method to setup a Python 2.7 environment. The Anaconda python distribution can be downloaded from [here](https://www.anaconda.com/products/distribution#Downloads) and installed in users' home directory. After Anaconda is installed, use these commands (one at a time) to generate, activate, and set up conda:

```bash
conda create -n flow python=2.7.17
```
```bash
conda activate flow
pip install cython six setuptools scons
```
Note the path anaconda is installed to, and update the paths for python in `buildoptions.py` :
```bash
PYTHON_INCDIR  = [usrhome+'/anaconda3/envs/flow/include/python2.7']
PYTHON_LIBSDIR = [usrhome+'/anaconda3/envs/flow/lib']
```
Delete `$HOME/anaconda3/envs/flow/lib/libgfortran.so` if it exists. (it is a symlink). It conflicts with the main `ligbfortran.so` library
## Metis
There is a precompiled libmetis included in the repo but it most likely will not work. ~~Download from [here](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz)~~. The original source link doesnt work anymore. The required version can be found via a search engine.
```bash
cd ~/temp
mkdir metis && cd $_
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
tar xf metis-4.0.3.tar.gz
cd metis-4.0.3
make
```
After building `libmetis.a` copy it to `external/x64-linux/metis-4.0`
## Verdandi
The source for verdandi is included in the repo. Verdandi depends on `liblua.a` . The source for liblua and a precompiled `liblua.a` file is also included inside the repo. Like `libmetis.a`, `liblua.a` will not be usable most likely and needs to be rebuilt. 
Finally, to build liblua go to `external/x64-linux/verdandi-1.7/include/lua` and run `make linux`. 

## Building CRIMSON flowsolver
After all the dependencies have been built, we can proceed to building the flowsolver. remember to activate the python environment. You need to update the paths to the dependencies in the `buildoptions.py` to point the locations where the dependencies have actually been installed. The paths that might need to set are identified with a `#setme` comment. Before running `scons`, activate the `flow` env that was created earlier by running `conda activate flow`.

Now, to build the code run `scons -j4` to use 4 cores for the build process. The build system accepts a few options:
1. debug: `scons debug=1` builds in debug mode, while `scons debug=0` builds the optimized version
2. test: `scons test=1` builds the test
3. profile: `scons profile=1` adds `-pg` to the build process to enable profiling using gprof.
All of these options can be specified at once.

# Running CRIMSON flowsolver
Before running the solver update LD_LIBRARY_PATH with the paths of the dependencies. Run `ldd bin/flowsolver | grep "not found"` to see if any libraries are not found. 

## Running tests
In the `testbin` folder there are example test runner scripts that can be copied and edited to run the tests. Copy/rename the testing script to `mytest` and update the variables `CRIMSON_FLOWSOLVER_HOME`, `PHASTA_CONFIG` and `LD_LIBRARY_PATH` in the script to match the present crimson installation and dependent libraries. Then begin the test by running

```bash
./mytest N 
```
where `N` is the number of process you want to use.