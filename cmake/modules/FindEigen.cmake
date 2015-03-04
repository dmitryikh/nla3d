# - Find Eigen
# Find the Eigen libraries
#
# This module defines the following variables:
#   EIGEN_FOUND -True if EIGEN_INCLUDE_DIR are found
#   EIGEN_INCLUDE_DIR - where to find Eigen.h.
#   EIGEN_INCLUDE_DIRS - set when EIGEN_INCLUDE_DIR found

include(FindPackageHandleStandardArgs)

set(EIGEN_ROOT "" CACHE PATH "Folder contains eigen")

# Find include dir
find_path(EIGEN_INCLUDE_DIR "Eigen/Core"
    PATHS ${EIGEN_ROOT}
    PATH_SUFFIXES  eigen3
    DOC "The path the the directory that contains Eigen.h")

find_package_handle_standard_args(Eigen DEFAULT_MSG
    EIGEN_INCLUDE_DIR)

if(EIGEN_FOUND)
    set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})
endif()
