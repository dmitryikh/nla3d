#   CMake module to find easylogging++ library
include(FindPackageHandleStandardArgs)
	
FIND_PATH(
    EASYLOGGINGPP_INCLUDE_DIR
    easylogging++.h
    NO_CMAKE_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/site-src/easyloggingpp/src
    DOC "Path to directory where easylogging++.h is stored"
)
	
FIND_PACKAGE_HANDLE_STANDARD_ARGS(EASYLOGGINGPP DEFAULT_MSG EASYLOGGINGPP_INCLUDE_DIR)
