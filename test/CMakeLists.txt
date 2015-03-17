# unit tests
add_executable(mat_math_test_exe mat_math.cpp ../src/sys.cpp)
target_link_libraries(mat_math_test_exe math)
add_test(NAME mat_math_test COMMAND mat_math_test_exe)



add_executable(mat_sparse_pardiso_test_exe mat_sparse_pardiso.cpp ../src/sys.cpp)
target_link_libraries(mat_sparse_pardiso_test_exe math ${MKL_LIBRARIES})
add_test(NAME mat_sparse_pardiso_test COMMAND mat_sparse_pardiso_test_exe)



add_executable(mat_prod_test_exe mat_prod/mat_prod_test.cpp ../src/sys.cpp)
if (nla3d_blas)
  target_link_libraries(mat_prod_test_exe math ${MKL_LIBRARIES})
else()
  target_link_libraries(mat_prod_test_exe math)
endif()
add_test(NAME mat_prod_test COMMAND mat_prod_test_exe -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)

add_executable(eigen_mat_prod_test_exe mat_prod/eigen_mat_prod_test.cpp ../src/sys.cpp)
if (nla3d_blas)
  target_link_libraries(eigen_mat_prod_test_exe math ${MKL_LIBRARIES})
else()
  target_link_libraries(eigen_mat_prod_test_exe math)
endif()
add_test(NAME eigen_mat_prod_test COMMAND eigen_mat_prod_test_exe -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)


add_test(NAME 3d_damper_test COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/3d_damper/model.cdb -material Neo-Hookean 10 5000 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/3d_damper/ansys/loading_curve_ansys.txt -reaction TOP_SIDE -threshold 3.0 -element SOLID81)
 

add_test(NAME a2000_damper_test COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/a2000_damper/a2000.cdb -element PLANE41 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/a2000_damper/ansys/loading_curve_ansys.txt -threshold 0.2 -reaction DISP_BC)

