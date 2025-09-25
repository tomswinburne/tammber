# FOR LAMMPS CMAKE BUILD

# set installation location
set(CMAKE_INSTALL_PREFIX "$ENV{PREFIX}" CACHE BOOL "" FORCE)
set(CMAKE_BUILD_PARALLEL_LEVEL 4) 
#set(Python_FIND_FRAMEWORK LAST)

#set(PYTHON_EXECUTABLE /path/to/python CACHE BOOL "" FORCE) 

# enforce c++11 standards
set(CCFLAGS -g -O3 -std=c++11)

# compile a binary and static library
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_LIB ON CACHE BOOL "" FORCE)

# allow error messages (very useful)
set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)

# minimal packages to run example (MANYBODY and ML-SNAP)
set(ALL_PACKAGES MANYBODY ML-SNAP)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
