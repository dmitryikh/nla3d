

add_subdirectory (math)

file (GLOB NLA3D_LIB_SOURCES "*.cpp" "materials/*.cpp" "elements/*.cpp")

add_library(nla3d_lib STATIC ${NLA3D_LIB_SOURCES})
target_link_libraries(nla3d_lib math ${MKL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set_target_properties(nla3d_lib PROPERTIES COTIRE_PREFIX_HEADER_INCLUDE_PATH
  "${CMAKE_SOURCE_DIR}/site-src")
cotire(nla3d_lib)
