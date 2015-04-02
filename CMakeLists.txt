cmake_minimum_required (VERSION 2.6)

project (NLA3D)

#options

set (nla3d_multithreaded OFF)
set (nla3d_debug OFF)
#NOTE: now SOLID80 can use blas, but results not converged.. need to fix
set (nla3d_blas OFF)
set(MKL_STATIC 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

if (nla3d_multithreaded)
  if(UNIX)
    find_package(Threads)
  endif()
  set(MKL_MULTI_THREADED ON)
  find_package(OpenMP)
  if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
else()
  set(MKL_MULTI_THREADED OFF)
endif() #nla3d_multithreaded

if (nla3d_debug)
  if(UNIX) 
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=gnu++11")
  endif()
  #todo make config=Debug by default for MSVC
else()
  if(UNIX) 
    # in case of release under linux use -O2 option
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O2 -msse2 -std=gnu++11")
  endif()
endif()

if (nla3d_blas)
  add_definitions( -DNLA3D_USE_BLAS)
endif()

if(WIN32)
  #this settings works for me.. But for you it needs to be changed
  set(MKL_ARCH "em64t")
else()
  set(MKL_ARCH "em64t")
endif()

find_package(MKL)
if (MKL_FOUND)
  include_directories(${MKL_INCLUDE_DIRS})
  link_directories(${MKL_LIBRARY_DIRS})
else()
  message(WARNING "Can't find MKL")
endif()


SET (CMAKE_COLOR_MAKEFILE ON)
enable_testing()

# site-src contains source codes of external libraries
SET (SITE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/site-src)

include_directories(${SITE_SOURCE_DIR}/eigen)
include_directories(${SITE_SOURCE_DIR}/easyloggingpp/src)
include_directories (src)

add_subdirectory (src/math)
add_subdirectory (test)


#TODO: make recursively add
file (GLOB NLA3D_LIB_SOURCES
	  "src/*.cpp" "src/materials/*.cpp" "src/elements/*cpp")

# here is a main target - nla3d.exe
set (NLA3D_APP_SOURCES "src/main.cpp")
list(REMOVE_ITEM NLA3D_LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp) 

#collect all sources except main.cpp into a static link
#this static link could be conviniently used for tests
add_library(nla3d_lib STATIC ${NLA3D_LIB_SOURCES})

target_link_libraries(nla3d_lib math ${MKL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})


#main executable of the project: nla3d.exe program
add_executable(nla3d ${NLA3D_APP_SOURCES})
target_link_libraries(nla3d nla3d_lib)


# Add logic to INSTALL.vcproj to copy app.exe to destination directory
install (TARGETS nla3d
		 RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin)

# taken from stackoverflow: How to copy DLL files into the same folder as the executable using CMake
add_custom_command(TARGET nla3d POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  ${MKL_IOMP5_RUNTIME_LIB}
  $<TARGET_FILE_DIR:nla3d>)

add_custom_command(TARGET nla3d POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  ${CMAKE_CURRENT_SOURCE_DIR}/files/logger.conf
  $<TARGET_FILE_DIR:nla3d>)
