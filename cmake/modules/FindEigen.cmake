# - Find Eigen
# Find the Eigen libraries

include(FindPackageHandleStandardArgs)
# Find include dir
find_path(EIGEN_INCLUDE_DIR "Eigen/Core"
    NO_CMAKE_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/site-src/eigen
    DOC "The path the the directory that contains Core, Eigen and others")

find_package_handle_standard_args(Eigen DEFAULT_MSG
    EIGEN_INCLUDE_DIR)
