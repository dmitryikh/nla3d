# Collect sources into the variable MATH_SOURCES with 
# having to explicitly list each header and source file
file (GLOB MATH_SOURCES
	  "*.cpp")

add_library(math STATIC ${MATH_SOURCES})

# Creates a folder "libraries" and adds target project (math.vcproj) under it
set_property(TARGET math PROPERTY FOLDER "libraries")

