cmake_minimum_required(VERSION 3.5)

set(PREFIX "$ENV{PREFIX}") # example value
# ENSURE ALL COMPILATIONS USE THE SAME MPI AS THE RUNTIME ENVIRONMENT!!
set(CMAKE_CXX_COMPILER mpicxx) # example value
set(CMAKE_C_COMPILER mpicc) # example value

# Setting VERBOSE=1 will generate a lot of output for debugging purposes
# To toggle LAMMPS output use the LogLammps option in input/ps-config.xml
set(VERBOSE 0)

# BOOST_LOG is often raises linking errors- disabled by default, log to stdout
set(USE_BOOST_LOG 0)

# If ON force search for static BOOST libs. Found to be essential for BlueGene
set(Boost_USE_STATIC_LIBS OFF)


set(CMAKE_PREFIX_PATH "${PREFIX}")
set(CMAKE_INSTALL_PREFIX "${PREFIX}")
set(CMAKE_INCLUDE_PATH ${PREFIX}/include)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib)


project (TAMMBER)

enable_testing()


# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

# Enable C++ and ensure C++11
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS ON) #...without compiler extensions like gnu++11

# Find python for testing
find_package(PythonInterp 2 REQUIRED)

# Boost: Identify components we need to load from Boost and load them
if(USE_BOOST_LOG)
	add_definitions(-DUSE_BOOST_LOG)
	set(BOOST_COMPONENTS random serialization timer system log log_setup filesystem thread regex chrono atomic date_time program_options)
else()
	set(BOOST_COMPONENTS random serialization timer system filesystem thread regex chrono atomic date_time program_options)
endif(USE_BOOST_LOG)


find_package(Boost 1.58.0 REQUIRED COMPONENTS ${BOOST_COMPONENTS})
foreach(mycomponent ${BOOST_COMPONENTS})
  list(APPEND BOOST_LIBRARIES "Boost::${mycomponent}")
endforeach(mycomponent)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Threads REQUIRED)

find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})

find_package(NAUTY REQUIRED)
include_directories(${NAUTY_INCLUDE_DIRS})

find_package(LAMMPS REQUIRED)
include_directories(${LAMMPS_INCLUDE_DIRS})

include_directories(core)
include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})

file(GLOB CORE_FILES "core/*")
file(GLOB TAMMBER_CORE_FILES "tammber-core/*")

set(EXTRA_FLAGS -g -DLAMMPS_EXCEPTIONS -DLAMMPS_LIB_MPI -DISOMORPHIC -DBOOST_BIND_GLOBAL_PLACEHOLDERS)
set(EXTRA_LINK_FLAGS pthread lammps_mpi)

set(TAMMBER_TARGETS tammber tammber-dbextract tammber-analyze tammber-neb tammber-md)
# must be subset of TAMMBER_TARGETS
set(FORCED_VERBOSE_TARGETS tammber-md tammber-neb)

add_executable(tammber ${CORE_FILES} ${TAMMBER_CORE_FILES} tammber/main.cpp)
add_executable(tammber-dbextract ${CORE_FILES} ${TAMMBER_CORE_FILES} tammber/dbextract.cpp)
add_executable(tammber-analyze ${CORE_FILES} ${TAMMBER_CORE_FILES} tammber/analyze.cpp)
add_executable(tammber-neb ${CORE_FILES} ${TAMMBER_CORE_FILES} tammber/neb.cpp)
add_executable(tammber-md ${CORE_FILES} ${TAMMBER_CORE_FILES} tammber/md.cpp)



if(VERBOSE)
	message("VERBOSE OUTPUT ON FOR ALL TARGETS")
	add_definitions(-DVERBOSE -g)
else()
	if(FORCED_VERBOSE_TARGETS)
		foreach(mytarget ${FORCED_VERBOSE_TARGETS})
			message("FORCING VERBOSITY FOR ${mytarget}")
			target_compile_definitions(${mytarget} PUBLIC VERBOSE)
		endforeach(mytarget)
	endif(FORCED_VERBOSE_TARGETS)
endif(VERBOSE)

if(EXTRA_FLAGS)
   message(STATUS "Compiling with these extra flags: ${EXTRA_FLAGS}")
   add_definitions(${EXTRA_FLAGS})
endif()

# Base dependencies for all targets
if(EXTRA_LINK_FLAGS)
   message(STATUS "Linking with these extra flags: ${EXTRA_LINK_FLAGS}")
endif()

set(MY_DEPENDENCIES ${CMAKE_THREAD_LIBS_INIT}
										${LAMMPS_LIBRARIES}
		    						${DB_LIBRARIES}
		    						${NAUTY_LIBRARIES}
		    						${BOOST_LIBRARIES}
		    						${MPI_CXX_LIBRARIES}
		    						${EXTRA_LINK_FLAGS}
)

# add dependencies for all targets
include(GNUInstallDirs)

foreach(mytarget ${TAMMBER_TARGETS})
  # add base dependencies
  target_include_directories( ${mytarget} PUBLIC tammber-core)
  target_link_libraries( ${mytarget} ${MY_DEPENDENCIES})
  install(TARGETS ${mytarget} DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach(mytarget)
