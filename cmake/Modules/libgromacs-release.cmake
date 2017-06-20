#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libgromacs" for configuration "Release"
set_property(TARGET libgromacs APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libgromacs PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/usr/lib/openmpi/lib/libmpi.so;/usr/lib/x86_64-linux-gnu/libhwloc.so;/usr/lib/x86_64-linux-gnu/libz.so;dl;rt;m;gmxfftw;/usr/lib/libblas.so.3gf;/usr/lib/liblapack.so.3gf;/usr/lib/libblas.so.3gf;-lpthread;-fopenmp"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgromacs_mpi.so.2.0.0"
  IMPORTED_SONAME_RELEASE "libgromacs_mpi.so.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS libgromacs )
list(APPEND _IMPORT_CHECK_FILES_FOR_libgromacs "${_IMPORT_PREFIX}/lib/libgromacs_mpi.so.2.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
