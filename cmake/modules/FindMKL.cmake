# - Find Intel MKL
# Find the MKL libraries
#

include(FindPackageHandleStandardArgs)

set(_INTEL_ROOT "/opt/intel" CACHE PATH "Folder contains intel libs")
set(_MKL_ROOT ${_INTEL_ROOT}/mkl CACHE PATH "Folder contains MKL")

# Find include dir
find_path(_MKL_INCLUDE_DIR mkl.h
    PATHS ${_MKL_ROOT}/include)


find_package_handle_standard_args(MKL DEFAULT_MSG
    _INTEL_ROOT _MKL_ROOT _MKL_INCLUDE_DIR)

if(MKL_FOUND)
    set(MKL_INCLUDE_DIR ${_MKL_INCLUDE_DIR})
    set(INTEL_ROOT ${_INTEL_ROOT})
    set(MKL_ROOT ${_MKL_ROOT})
endif()
