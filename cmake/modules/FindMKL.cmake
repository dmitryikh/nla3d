# - Find Intel MKL
# Find the MKL libraries
#
# Options:
#
#   MKL_STATIC       :   use static linking
#   MKL_MULTI_THREADED:   use multi-threading
#   MKL_SDL           :   Single Dynamic Library interface
#   MKL_ARCH          : binary architecture {32, 64, em64t, ... }
#
# This module defines the following variables:
#
#   MKL_FOUND            : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDE_DIR      : where to find mkl.h, etc.
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.
#   MKL_IOMP5_RUNTIME_LIB


include(FindPackageHandleStandardArgs)

set(INTEL_ROOT "/opt/intel" CACHE PATH "Folder contains intel libs")
set(MKL_ROOT ${INTEL_ROOT}/mkl CACHE PATH "Folder contains MKL")

# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h
    PATHS ${MKL_ROOT}/include)


# Find libraries
# Handle suffix
set(_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(WIN32)
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib)
    endif()
else()
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
    endif()
endif()


# MKL is composed by four layers: Interface, Threading, Computational and RTL
if (LINUX)
  set(MKL_LIB_DIR ${MKL_ROOT}/lib/${MKL_ARCH}/)
elseif(APPLE)
  set(MKL_LIB_DIR ${MKL_ROOT}/lib/)
else()
  set(MKL_LIB_DIR ${MKL_ROOT}/${MKL_ARCH}/lib)
endif()

if(MKL_SDL)
    find_library(MKL_LIBRARY mkl_rt
        PATHS ${MKL_LIB_DIR})
    set(MKL_MINIMAL_LIBRARY ${MKL_LIBRARY})
else()
    ######################### Interface layer #######################
    if(WIN32)
      find_library(MKL_INTERFACE_LIBRARY mkl_intel_lp64 PATHS ${MKL_LIB_DIR})
      find_library(MKL_INTERFACE_LIBRARY mkl_intel_c PATHS ${MKL_LIB_DIR})
    else()
      find_library(MKL_INTERFACE_LIBRARY mkl_intel_lp64 PATHS ${MKL_LIB_DIR})
      find_library(MKL_INTERFACE_LIBRARY mkl_intel PATHS ${MKL_LIB_DIR})
    endif()


    ######################## Threading layer ########################
    if(MKL_MULTI_THREADED)
      if (CMAKE_COMPILER_IS_GNUCXX)
        set(MKL_THREADING_LIBNAME mkl_gnu_thread)
      else()
        set(MKL_THREADING_LIBNAME mkl_intel_thread)
      endif()
    else()
        set(MKL_THREADING_LIBNAME mkl_sequential)
    endif()

    find_library(MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME}
      PATHS ${MKL_LIB_DIR})

    ####################### Computational layer #####################
    find_library(MKL_CORE_LIBRARY mkl_core
      PATHS ${MKL_LIB_DIR})

    ############################ RTL layer ##########################
    if(WIN32)
        set(MKL_RTL_LIBNAME libiomp5md)
    else()
        set(MKL_RTL_LIBNAME iomp5)
    endif()
    find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
        PATHS ${INTEL_ROOT}/lib)
    if(CMAKE_COMPILER_IS_GNUCXX)
      #need to group undex g++ in urder to resolve all references
      set(MKL_LIBRARY "-Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY}${MKL_SCALAPACK_LIBRARY} ${MKL_RTL_LIBRARY} ${MKL_SOLVER_LIBRARY} -Wl,--end-group")
    else()
      set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_RTL_LIBRARY} ${MKL_SOLVER_LIBRARY})
    endif()
    set(MKL_MINIMAL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
endif()


if(WIN32)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .dll)
  find_library(MKL_IOMP5_RUNTIME_LIB libiomp5md
      PATHS ${MKL_ROOT}/${MKL_ARCH}/bin/)
elseif(APPLE)
  set(MKL_RTL_LIBNAME iomp5)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib)
  find_library(MKL_IOMP5_RUNTIME_LIB ${MKL_RTL_LIBNAME}
      PATHS ${INTEL_ROOT}/lib)
else()
  set(MKL_RTL_LIBNAME iomp5)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
  find_library(MKL_IOMP5_RUNTIME_LIB ${MKL_RTL_LIBNAME}
      PATHS ${MKL_LIB_DIR})
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})


find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_INCLUDE_DIR MKL_LIBRARY MKL_MINIMAL_LIBRARY MKL_IOMP5_RUNTIME_LIB)

if(MKL_FOUND)
    set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
    set(MKL_LIBRARIES ${MKL_LIBRARY})
    set(MKL_MINIMAL_LIBRARIES ${MKL_LIBRARY})
endif()
