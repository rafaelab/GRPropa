cmake_minimum_required(VERSION 2.6)
project(GRPropa C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(GRPROPA_EXTRA_SOURCES)
set(GRPROPA_EXTRA_INCLUDES)
set(GRPROPA_EXTRA_LIBRARIES)
set(GRPROPA_SWIG_DEFINES)

if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
	message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

# Set default build-type to release to enable performance improvements 
if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
MESSAGE(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")


# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# kiss (provided)
add_subdirectory(libs/kiss)
list(APPEND GRPROPA_EXTRA_LIBRARIES kiss)
list(APPEND GRPROPA_EXTRA_INCLUDES libs/kiss/include)

# HepID (provided)
add_subdirectory(libs/HepPID)
list(APPEND GRPROPA_EXTRA_LIBRARIES HepPID)
list(APPEND GRPROPA_EXTRA_INCLUDES libs/HepPID/include)

# SAGA
option(ENABLE_SAGA "Compiling with SAGA" OFF)
if(ENABLE_SAGA)
	find_package(SAGA)
	if(SAGA_FOUND)
		list(APPEND GRPROPA_EXTRA_INCLUDES ${SAGA_INCLUDE_DIR})
		list(APPEND GRPROPA_EXTRA_LIBRARIES ${SAGA_LIBRARY})
		add_definitions (-DGRPROPA_HAVE_SAGA)
		list(APPEND GRPROPA_SWIG_DEFINES -DGRPROPA_HAVE_SAGA)
		list(APPEND GRPROPA_SWIG_DEFINES -I${SAGA_INCLUDE_DIR})
	endif(SAGA_FOUND)

	# SQLite 3 (needed for SAGA)
	find_path(SQLITE3_INCLUDE_DIR sqlite3.h)
	find_library(SQLITE3_LIBRARY libsqlite3)
	include_directories(${SQLITE3_INCLUDE_DIR})
	set(SQLITE3_CFLAGS "")
	if(SQLITE3_LIBRARY AND SQLITE3_INCLUDE_DIR)
	    message(STATUS "SQLite3 found!")
	    message(STATUS "  Include:        ${SQLITE3_INCLUDE_DIR}")
	    message(STATUS "  Library:        ${SQLITE3_LIBRARY}")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SQLITE3_CFLAGS}")
	    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ${SQLITE3_LIBRARY} ")
	    list(APPEND GRPROPA_EXTRA_INCLUDES ${SQLITE3_INCLUDE_DIR})
	    list(APPEND GRPROPA_EXTRA_LIBRARIES ${SQLITE3_LIBRARY})
	else()
	    message(STATUS "SQLite3: NOT Found!")    
	endif()
endif(ENABLE_SAGA)

# OpenMP (optional for shared memory multiprocessing)
option(ENABLE_OPENMP "OpenMP for multithreading" ON)
if(ENABLE_OPENMP)
include(FindOpenMP)
	if(OPENMP_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif(OPENMP_FOUND)
endif(ENABLE_OPENMP)

# FFTW3F (optional for turbulent magnetic fields)
find_package(FFTW3F)
if(FFTW3F_FOUND)
	list(APPEND GRPROPA_EXTRA_INCLUDES ${FFTW3F_INCLUDE_DIR})
	list(APPEND GRPROPA_EXTRA_LIBRARIES ${FFTW3F_LIBRARY})
	add_definitions(-DGRPROPA_HAVE_FFTW3F)
	list(APPEND GRPROPA_SWIG_DEFINES -DGRPROPA_HAVE_FFTW3F)
endif(FFTW3F_FOUND)



# ----------------------------------------------------------------------------
# Library and Binary
# ----------------------------------------------------------------------------
include_directories(include ${GRPROPA_EXTRA_INCLUDES})

add_library(grpropa SHARED
	src/Random.cpp
	src/Clock.cpp
	src/ModuleList.cpp
	src/Module.cpp
	src/Candidate.cpp
	src/ParticleState.cpp
	src/ProgressBar.cpp
	src/Cosmology.cpp
	src/Source.cpp
	src/Common.cpp
	src/PhotonBackground.cpp
	src/GridTools.cpp
	src/module/BreakCondition.cpp
	src/module/Boundary.cpp
	src/module/Observer.cpp
	src/module/SimplePropagation.cpp
	src/module/PropagationCK.cpp
	src/module/InverseCompton.cpp
	src/module/PairProduction.cpp
    src/module/Synchrotron.cpp
	src/module/Redshift.cpp
	src/module/Output.cpp
	src/module/TextOutput.cpp
	src/module/Tools.cpp
	src/magneticField/MagneticField.cpp
	src/magneticField/MagneticFieldGrid.cpp
	src/magneticField/TurbulentMagneticField.cpp
	src/magneticField/JF12Field.cpp
	${GRPROPA_EXTRA_SOURCES}
)
target_link_libraries(grpropa ${GRPROPA_EXTRA_LIBRARIES})


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
add_definitions(-DGRPROPA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS grpropa DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY data/ DESTINATION share/grpropa/ PATTERN ".git" EXCLUDE)

# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via SWIG" ON)
if(ENABLE_PYTHON)
	FIND_PACKAGE(SWIG REQUIRED)
	include(python/Python.cmake)
	include_directories(${PYTHON_INCLUDE_PATH})

	file(GLOB_RECURSE GRPROPA_INCLUDES include/*.h)
	set_source_files_properties( ${CMAKE_CURRENT_BINARY_DIR}/grpropa_wrap.cxx PROPERTIES GENERATED true )
	add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/grpropa_wrap.cxx
		COMMAND swig -c++ -python -I${CMAKE_SOURCE_DIR}/include ${GRPROPA_SWIG_DEFINES} -o ${CMAKE_CURRENT_BINARY_DIR}/grpropa_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/python/grpropa.i
		DEPENDS ${CMAKE_SOURCE_DIR}/python/grpropa.i ${GRPROPA_INCLUDES} )

	add_library(grpropa-swig MODULE ${CMAKE_CURRENT_BINARY_DIR}/grpropa_wrap.cxx)
	set_target_properties(grpropa-swig PROPERTIES PREFIX "")
	set_target_properties(grpropa-swig PROPERTIES OUTPUT_NAME "_grpropa")
	target_link_libraries(grpropa-swig grpropa ${PYTHON_LIBRARIES})

	install(DIRECTORY "${CMAKE_SOURCE_DIR}/python/grpropa" DESTINATION ${PYTHON_SITE_PACKAGES})
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/grpropa.py" DESTINATION ${PYTHON_SITE_PACKAGES}/grpropa)
	install(TARGETS grpropa-swig LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}/grpropa)
endif(ENABLE_PYTHON)
