cmake_minimum_required (VERSION 2.8.10)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

########################### Tenstream Project: ###################################
project (maxrand)
enable_language(C Fortran)


# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS "-cpp -ffree-line-length-none -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbacktrace -finit-real=nan -W -Wall -Wuninitialized -g -pg -fcheck=all -fbounds-check -pedantic -Wsurprising -ffpe-trap=invalid,zero,overflow -fno-range-check")

elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS "-cpp -traceback -extend_source -g -sox -no-wrap-margin ")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model source")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")
else (Fortran_COMPILER_NAME MATCHES "*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Set Output Directories
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp archives")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp libraries")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Output directory for Dohp executables")

# Add Lapack
find_package( LAPACK )
if(${LAPACK_FOUND})
  list (APPEND EXTLIB "${LAPACK_LIBRARIES}")
endif(${LAPACK_FOUND})

# Add Source Directories
add_subdirectory (src)

find_program(CTAGS_PATH ctags)
if(CTAGS_PATH)
    message(STATUS "Found ctags: ${CTAGS_PATH}")
    set_source_files_properties( tags PROPERTIES GENERATED true)
    add_custom_target(tags ALL
        COMMAND ctags -R ${CMAKE_CURRENT_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Creating ctags" VERBATIM
        )
endif(CTAGS_PATH)
