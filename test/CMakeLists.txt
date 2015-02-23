# unit tests
add_executable(mat_math_test_exe mat_math.cpp ../src/sys.cpp)
target_link_libraries(mat_math_test_exe math)
add_test(NAME mat_math_test COMMAND mat_math_test_exe)



add_executable(mkl_band_test_exe mkl_band.cpp ../src/sys.cpp)
target_link_libraries(mkl_band_test_exe math mkl_intel_c mkl_intel_thread libiomp5md mkl_core mkl_solver)
add_test(NAME mkl_band_test COMMAND mkl_band_test_exe)



add_executable(mat_sparse_pardiso_test_exe mat_sparse_pardiso.cpp ../src/sys.cpp)
target_link_libraries(mat_sparse_pardiso_test_exe math mkl_intel_c mkl_intel_thread libiomp5md mkl_core mkl_solver)
add_test(NAME mat_sparse_pardiso_test COMMAND mat_sparse_pardiso_test_exe)



add_executable(mat_prod_test_exe mat_prod/mat_prod_test.cpp ../src/sys.cpp)
target_link_libraries(mat_prod_test_exe math)
add_test(NAME mat_prod_test COMMAND mat_prod_test_exe -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)



# test cases
# build exe for test 3d_damper
add_executable(3d_damper_test_exe 3d_damper/3d_damper_test.cpp)
target_link_libraries(3d_damper_test_exe nla3d_lib)
add_test(NAME 3d_damper_test COMMAND 3d_damper_test_exe ${PROJECT_SOURCE_DIR}/test/3d_damper/model.cdb -material Neo-Hookean 10 5000 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/3d_damper/ansys/loading_curve_ansys.txt)
# taken from stackoverflow: How to copy DLL files into the same folder as the executable using CMake
add_custom_command(TARGET 3d_damper_test_exe POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  ${MKL_PATH}/ia32/bin/libiomp5md.dll
  $<TARGET_FILE_DIR:3d_damper_test_exe>)
 

# build exe for test 3d_damper
add_executable(a2000_damper_test_exe a2000_damper/a2000_damper_test.cpp)
target_link_libraries(a2000_damper_test_exe nla3d_lib)
add_test(NAME a2000_damper_test COMMAND a2000_damper_test_exe ${PROJECT_SOURCE_DIR}/test/a2000_damper/a2000.cdb -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/a2000_damper/ansys/loading_curve_ansys.txt)
# taken from stackoverflow: How to copy DLL files into the same folder as the executable using CMake
add_custom_command(TARGET a2000_damper_test_exe POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  ${MKL_PATH}/ia32/bin/libiomp5md.dll
  $<TARGET_FILE_DIR:a2000_damper_test_exe>)

