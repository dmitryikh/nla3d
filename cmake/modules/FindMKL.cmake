# - Find Intel MKL
# Find the MKL libraries
#

include(FindPackageHandleStandardArgs)

set(_INTEL_ROOT "/opt/intel" CACHE PATH "Folder contains intel libs")
set(_MKL_ROOT ${_INTEL_ROOT}/mkl CACHE PATH "Folder contains MKL")

# Find include dir
find_path(_MKL_INCLUDE_DIR mkl.h
    PATHS ${_MKL_ROOT}/include)

# Find library dir
find_path(_MKL_LIBS_DIR
    NAMES libmkl_core.a mkl_core.lib
    PATHS ${_MKL_ROOT}/lib
          ${_MKL_ROOT}/lib/intel64)

find_package_handle_standard_args(MKL DEFAULT_MSG
    _INTEL_ROOT _MKL_ROOT _MKL_INCLUDE_DIR _MKL_LIBS_DIR)


if(MKL_FOUND)
    set(MKL_INCLUDE_DIR ${_MKL_INCLUDE_DIR})
    set(MKL_LIBS_DIR ${_MKL_LIBS_DIR})
    set(INTEL_ROOT ${_INTEL_ROOT})
    set(MKL_ROOT ${_MKL_ROOT})
    set(MKL_LINKER_FLAGS " -L${MKL_LIBS_DIR} -Wl,-rpath,${MKL_LIBS_DIR}")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(MKL_LINKER_FLAGS "${MKL_LINKER_FLAGS},--no-as-needed")
    endif()
endif()
