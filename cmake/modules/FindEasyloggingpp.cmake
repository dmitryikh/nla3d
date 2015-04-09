#   CMake module to find easylogging++ library
	
FIND_PATH(
    EASYLOGGINGPP_INCLUDE_DIRS
    NAMES src/easylogging++.h
    HINTS $ENV{EASYLOGGINGPP_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/site-src
    PATHS /usr/local/include
          /usr/include
)
	
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(EASYLOGGINGPP DEFAULT_MSG EASYLOGGINGPP_INCLUDE_DIRS)
