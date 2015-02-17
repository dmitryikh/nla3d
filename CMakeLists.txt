cmake_minimum_required (VERSION 2.6)

# Maps to a solution file (Tutorial.sln). The solution will 
# have all targets (exe, lib, dll) as projects (.vcproj)
project (NLA3D)

# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
#set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Set compiler flags and options. 
# Here it is setting the Visual Studio warning level to 4
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")

# Command to output information to the console
# Useful for displaying errors, warnings, and debugging
message ("cxx Flags: " ${CMAKE_CXX_FLAGS})

# Sub-directories where more CMakeLists.txt exist
add_subdirectory (src/math)

# Properties->C/C++->General->Additional Include Directories
#include_directories ("${PROJECT_SOURCE_DIR}/math")
include_directories (src)


#now it's hardcoded because I don't know how to do it more flexible
set(MKL_PATH "C:/Documents and Settings/admin/My Documents/10.0.2.019/10.0.2.019")
include_directories (${MKL_PATH}/include)

#MKL libraries that needed to link with nla3d
add_library(mkl_intel_c STATIC IMPORTED)
set_property(TARGET mkl_intel_c PROPERTY IMPORTED_LOCATION ${MKL_PATH}/ia32/lib/mkl_intel_c.lib)

add_library(mkl_intel_thread STATIC IMPORTED)
set_property(TARGET mkl_intel_thread PROPERTY IMPORTED_LOCATION ${MKL_PATH}/ia32/lib/mkl_intel_thread.lib)

add_library(libiomp5md STATIC IMPORTED)
set_property(TARGET libiomp5md PROPERTY IMPORTED_LOCATION ${MKL_PATH}/ia32/lib/libiomp5md.lib)

add_library(mkl_core STATIC IMPORTED)
set_property(TARGET mkl_core PROPERTY IMPORTED_LOCATION ${MKL_PATH}/ia32/lib/mkl_core.lib)

add_library(mkl_solver STATIC IMPORTED)
set_property(TARGET mkl_solver PROPERTY IMPORTED_LOCATION ${MKL_PATH}/ia32/lib/mkl_solver.lib)

file (GLOB APP_SOURCES
	  "src/*.cpp")


# Set Properties->General->Configuration Type to Application(.exe)
# Creates nla3d.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer
add_executable (nla3d ${APP_SOURCES})

# Properties->Linker->Input->Additional Dependencies
target_link_libraries(nla3d math mkl_intel_c mkl_intel_thread libiomp5md mkl_core mkl_solver)

# Creates a folder "executables" and adds target 
# project (app.vcproj) under it
#set_property(TARGET app PROPERTY FOLDER "executables")

# Adds logic to INSTALL.vcproj to copy app.exe to destination directory
install (TARGETS nla3d
		 RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin)

# taken from stackoverflow: How to copy DLL files into the same folder as the executable using CMake
add_custom_command(TARGET nla3d POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  ${MKL_PATH}/ia32/bin/libiomp5md.dll
  $<TARGET_FILE_DIR:nla3d>)
