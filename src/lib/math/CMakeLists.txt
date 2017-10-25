# Collect sources into the variable MATH_SOURCES with 
# having to explicitly list each header and source file
file (GLOB MATH_SOURCES "*.cpp")

add_library(math STATIC ${MATH_SOURCES} ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)

# Creates a folder "libraries" and adds target project (math.vcproj) under it
set_property(TARGET math PROPERTY FOLDER "libraries")
set_target_properties(math PROPERTIES COTIRE_PREFIX_HEADER_INCLUDE_PATH
  "${CMAKE_SOURCE_DIR}/site-src")
cotire(math)

