cmake_minimum_required (VERSION 2.6)

project (NLA3D)
# Use mine findFoo.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

#options
set (nla3d_use_MKL ON)
set (nla3d_multithreaded OFF)
set (nla3d_debug OFF)
#TODO: now SOLID81 can use blas, but results not converged.. need to fix
# Do not use nla3d_blas for now..
set (nla3d_blas OFF)


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

if (nla3d_use_MKL)
    add_definitions( -DNLA3D_USE_MKL)
    # we can use blas procedures only if MKL library is used
    if (nla3d_blas)
      add_definitions( -DNLA3D_USE_BLAS)
    endif()

    find_package(MKL)
    if (MKL_FOUND)
      include_directories(${MKL_INCLUDE_DIR})

      # below code works fine only on UNIX systems..
      if (MKL_MULTI_THREADED)
        SET(MKL_LINKER_ARGS " -L${MKL_ROOT}/lib -Wl,-rpath,${MKL_ROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")
      else()
        SET(MKL_LINKER_ARGS " -L${MKL_ROOT}/lib -Wl,-rpath,${MKL_ROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
      endif()
        SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${MKL_LINKER_ARGS}" )

    else()
      message(WARNING "Can't find MKL")
    endif()
endif() # nla3d_use_MKL


find_package(EASYLOGGINGPP)
if (EASYLOGGINGPP_FOUND)
  include_directories(${EASYLOGGINGPP_INCLUDE_DIR})
else()
  message(WARNING "Can't find Easylogging++")
endif()

find_package(EIGEN)
if (EIGEN_FOUND)
  include_directories(${EIGEN_INCLUDE_DIR})
else()
  message(WARNING "Can't find Eigen")
endif()


SET (CMAKE_COLOR_MAKEFILE ON)
enable_testing()

include_directories (src)
include_directories (src/lib)

add_subdirectory (src/lib)
add_subdirectory (test)


# here is a main target - nla3d.exe
set (NLA3D_APP_SOURCES "src/main.cpp")

# an executable to use TRUSS3 element
set (NLA3D_TRUSS_APP_SOURCES "src/main_truss.cpp")

#main executable of the project: nla3d.exe program
add_executable(nla3d ${NLA3D_APP_SOURCES})
target_link_libraries(nla3d nla3d_lib)


add_executable(nla3d_truss ${NLA3D_TRUSS_APP_SOURCES})
target_link_libraries(nla3d_truss nla3d_lib)

# Add logic to INSTALL.vcproj to copy app.exe to destination directory
install (TARGETS nla3d
		 RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin)

# taken from stackoverflow: How to copy DLL files into the same folder as the executable using CMake
# add_custom_command(TARGET nla3d POST_BUILD
#  COMMAND ${CMAKE_COMMAND} -E copy_if_different
#  ${MKL_IOMP5_RUNTIME_LIB}
#  $<TARGET_FILE_DIR:nla3d>)
