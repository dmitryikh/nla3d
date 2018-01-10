cmake_minimum_required (VERSION 2.6)

project (NLA3D)

include(GNUInstallDirs)

# Use mine findFoo.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# cotire module will speed up compilation time
# see https://github.com/sakra/cotire
include(cotire)

#options
set (NLA3D_USE_MKL OFF
     CACHE BOOL "Use MKL library for math routines")
set (NLA3D_PYTHON OFF
     CACHE BOOL "Build python bindings")
set (nla3d_multithreaded OFF
    CACHE BOOL "Compile with OpenMP (not supported now)" FORCE)
#TODO: now SOLID81 can use blas, but results not converged.. need to fix
# Do not use NLA3D_BLAS for now..
set (NLA3D_BLAS OFF
     CACHE BOOL "Use blas routines for matrix manipulations (not worked properly)")

# TODO:
# CMAKE_BUILD_TYPE empty or Release means the same for us
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse2")

# TODO: this should depend on compiler
if(UNIX) 
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -fPIC")
endif()

# Linking on Mac Os was not working without this, see `cmake  --help-policy CMP0042`
if (APPLE)
       set(CMAKE_MACOSX_RPATH ON)
endif (APPLE)

if (nla3d_multithreaded)
  if(UNIX)
    find_package(Threads)
  endif()
  set(MKL_MULTI_THREADED ON)
  find_package(OpenMP)
  if (OPENMP_FOUND)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
else()
  set(MKL_MULTI_THREADED OFF)
endif() #nla3d_multithreaded


if (NLA3D_USE_MKL)
    add_definitions( -DNLA3D_USE_MKL)
    # we can use blas procedures only if MKL library is used
    if (NLA3D_BLAS)
      add_definitions( -DNLA3D_USE_BLAS)
    endif()

    find_package(MKL)
    if (MKL_FOUND)
      include_directories(${MKL_INCLUDE_DIR})

      # below code works fine only on UNIX systems..
      if (MKL_MULTI_THREADED)
          SET(MKL_LINKER_ARGS "${MKL_LINKER_FLAGS} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")
      else()
        SET(MKL_LINKER_ARGS "${MKL_LINKER_FLAGS} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
      endif()
      SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${MKL_LINKER_ARGS}" )
      SET( CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} ${MKL_LINKER_ARGS}" )

    else()
      message(WARNING "Can't find MKL")
    endif()
endif() # NLA3D_USE_MKL


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
# now we can use `make check` to build and run tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -L FUNC --output-on-failure)

include_directories (src)
include_directories (src/lib)

add_subdirectory (src/lib)


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

# install headers to install path
install(DIRECTORY "src/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nla3d"
          FILES_MATCHING PATTERN "*.h")

# generate python wrapper by SWIG. `swig` rule to generate wrapper sources, `nla3d_py` - to build library
IF (NLA3D_PYTHON)
    set(Python_ADDITIONAL_VERSIONS 2.7)
    FIND_PACKAGE(PythonInterp)
    FIND_PACKAGE(PythonLibs)
    FIND_PROGRAM(SWIG_EXECUTABLE NAMES swig3.0 swig)
    IF(SWIG_EXECUTABLE)
        MESSAGE("-- Found swig on \"${SWIG_EXECUTABLE}\"")
    ELSE()
        MESSAGE(WARNING "swig isn't found, you can't build nla3d_py")
    ENDIF()
    IF(PYTHONLIBS_FOUND AND SWIG_EXECUTABLE)
        MAKE_DIRECTORY(${CMAKE_BINARY_DIR}/python)
        ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_BINARY_DIR}/python/nla3d_wrapper.cxx
                           COMMAND ${SWIG_EXECUTABLE} -python -c++ -I${CMAKE_SOURCE_DIR}/src
                           -I${CMAKE_SOURCE_DIR}/src/lib -outdir ${CMAKE_BINARY_DIR}/python
                           -o python/nla3d_wrapper.cxx ${CMAKE_SOURCE_DIR}/python/nla3d.swg
                           DEPENDS ${CMAKE_SOURCE_DIR}/python/nla3d.swg)

        ADD_CUSTOM_TARGET(swig DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/python/nla3d_wrapper.cxx)

        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
        LINK_DIRECTORIES(${PYTHON_LIBRARIES})

        ADD_LIBRARY(nla3d_py SHARED EXCLUDE_FROM_ALL "${CMAKE_BINARY_DIR}/python/nla3d_wrapper.cxx")
        TARGET_LINK_LIBRARIES(nla3d_py ${PYTHON_LIBRARIES} nla3d_static_lib math)
        SET_TARGET_PROPERTIES(nla3d_py PROPERTIES OUTPUT_NAME "_nla3d"
                                                  PREFIX ""
                                                  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")
        if ( MSVC )
            SET_TARGET_PROPERTIES( nla3d_py PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
                                                       RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/python"
                                                       RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/python"
                                                       SUFFIX ".pyd")
        # etc for the other available configuration types (MinSizeRel, RelWithDebInfo)
        endif ( MSVC )
        if (APPLE)
            SET_TARGET_PROPERTIES(nla3d_py PROPERTIES SUFFIX ".so")
        endif (APPLE)
        ADD_DEPENDENCIES(nla3d_py swig)

        # run python tests with `checkpy`
        ADD_CUSTOM_TARGET(checkpy COMMAND ${CMAKE_CTEST_COMMAND} -L PYTHON --output-on-failure)
        ADD_DEPENDENCIES(checkpy nla3d_py)
    ENDIF()
ENDIF()

add_subdirectory (test EXCLUDE_FROM_ALL)
# taken from stackoverflow: How to copy DLL files into the same folder as the executable using CMake
# add_custom_command(TARGET nla3d POST_BUILD
#  COMMAND ${CMAKE_COMMAND} -E copy_if_different
#  ${MKL_IOMP5_RUNTIME_LIB}
#  $<TARGET_FILE_DIR:nla3d>)

