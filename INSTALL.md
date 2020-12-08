# Building TAMMBER / ParSplice

The following closely follows the ParSplice [wiki](https://gitlab.com/exaalt/parsplice/-/wikis/building-parsplice),

## Requirements

- Compiler with full C++11 support
- cmake 3.5 or higher

**Take care to ensure to use the same compiler is used throughout**

This is particularly relevant when using precompiled libraries (e.g. `module load`)


## Dependencies

### Prepare environment

- Set `PREFIX`, where dependencies will be installed, e.g. `${HOME}/.local`
	```bash
	export PREFIX=${HOME}/.local
	```

- Make destination directory
	```bash
	mkdir ${PREFIX}
	mkdir ${PREFIX}/lib
	mkdir ${PREFIX}/include
	```

- For testing on local machines it can be possible to install dependencies apart from `LAMMPS` with `apt`

	```bash
	sudo apt install libboost-all-dev libeigen3-dev libnauty-dev
	```
	In general, follow the instructions below.


### Boost:

- **Boost compilation is the most common source of problems, but is often available precompiled**

- Repeat: **Take care to ensure to use the same compiler is used throughout**

- Have a look with e.g. `module avail boost` we require a version `>=1.58`

- If available, install with e.g. `module load boost/1.58` and move on...

- If not available, unpack the archive:

	```bash
	cd deps
	tar -xvf boost_1_64_0.tar.gz
	cd boost*
	```

- Bootstrap to set your install prefix:

	```bash
	./bootstrap.sh --prefix=${PREFIX}
	```

- Build with e.g. gcc and c++11 support:

	```bash
	./b2 toolset=gcc cxxflags=-std=c++11 # gcc is example value
	```

- Install:

	```bash
	./b2 install
	```


### Eigen

- Have a look with e.g. `module avail eigen` we require a version `>=3`

- If available, install with e.g. `module load eigen/3.0.0` and move on...

- If not, unpack the archive:

	```bash
	cd deps
	tar -xvf eigen-eigen-*.tar.bz2
	cd eigen*
	```

- Build with `cmake` and install:

	```bash
	mkdir build
	cd build
	FC="nofortran" cmake -DCMAKE_INSTALL_PREFIX=${PREFIX}  ../
	make install
	```

	If you have a fortran compiler installed, `FC="nofortran"` will not be required.


### Nauty

- Have a look with e.g. `module avail nauty`

- If available, install with e.g. `module load nauty` and move on...

- If not, unpack the archive:

	```bash
	cd deps
	tar -xvf nauty*.tar.bz2
	cd nauty*
	```

- Configure, make and install:

	```bash
	./configure --prefix=${PREFIX}
	make
	cp nauty.a ${PREFIX}/lib/libnauty.a
	mkdir ${PREFIX}/include/nauty
	cp *.h ${PREFIX}/include/nauty
	```


### LAMMPS

- Download the `LAMMPS` repository:
	`git clone https://github.com/lammps/lammps.git`

- As of 01/09/2020 `lammps_gather()` and `lammps_scatter()` are now part of the `LAMMPS` library interface. You must therefore use the most recent version, obtained with e.g
```bash
git clone https://github.com/lammps/lammps.git
```
or replace `src/library.*` in your `LAMMPS` distribution

- Configure the `LAMMPS` build :
	```bash
	cd /path/to/lammps/src
	make yes-manybody # EAM potentials etc
	make yes-..... # whatever you want
	```
- Edit the makefile (`src/MAKE/Makefile.mpi` in this example) :
	```make
	...
	#CCFLAGS =	-g -O3 # previous
	CCFLAGS =	-g -O3 -std=c++11 # our modification
	...
	# LMP_INC =	-DLAMMPS_GZIP -DLAMMPS_MEMALIGN=64 # previous
	LMP_INC =	-DLAMMPS_GZIP -DLAMMPS_MEMALIGN=64  -DLAMMPS_EXCEPTIONS # our modification
	```

- Compile as a library and install :
	```bash
	make -j4 mpi mode=lib
	cp liblammps_mpi.a ${PREFIX}/lib
	mkdir ${PREFIX}/include/lammps
	cp *.h  ${PREFIX}/include/lammps/
	```
	where `-j4` will parallelize over 4 processors/threads

## TAMMBER
- Build binaries with `cmake` then `make`:
	```bash
	mkdir build
	cd build
	echo ${PREFIX} # make sure PREFIX is in the environment!
	cmake ../
	make -j4
	```
	where `-j4` will parallelize over 4 processors/threads. If this fails, check the values of `CMAKE_CXX_COMPILER` and `CMAKE_C_COMPILER` to ensure they correspond to the same as that used to compile the dependencies. Good luck!
