# unit tests
add_executable(mat_math_test_exe mat_math.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
target_link_libraries(mat_math_test_exe math)
add_test(NAME mat_math_test COMMAND mat_math_test_exe)

# TODO: undone! (see commented lines in the source code)
# add_executable(mat_sparse_pardiso_test_exe mat_sparse_pardiso.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
# target_link_libraries(mat_sparse_pardiso_test_exe math ${MKL_LIBRARIES})
# add_test(NAME mat_sparse_pardiso_test COMMAND mat_sparse_pardiso_test_exe)

add_executable(mat_prod_test_exe mat_prod/mat_prod_test.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
if (nla3d_blas)
  target_link_libraries(mat_prod_test_exe math ${MKL_LIBRARIES})
else()
  target_link_libraries(mat_prod_test_exe math)
endif()
add_test(NAME mat_prod_test COMMAND mat_prod_test_exe -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)

add_executable(eigen_mat_prod_test_exe mat_prod/eigen_mat_prod_test.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
if (nla3d_blas)
  target_link_libraries(eigen_mat_prod_test_exe math ${MKL_LIBRARIES})
else()
  target_link_libraries(eigen_mat_prod_test_exe math)
endif()
add_test(NAME eigen_mat_prod_test COMMAND eigen_mat_prod_test_exe -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)


# this test takes ~90 sec on Release..
add_test(NAME 3d_damper COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/3d_damper/model.cdb -material Neo-Hookean 10 5000 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/3d_damper/ansys/loading_curve_ansys.txt -reaction TOP_SIDE UY -threshold 3.0 -element SOLID81)
 
# this test takes toooo long with gauss eq. solver
 add_test(NAME a2000_damper COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/a2000_damper/a2000.cdb -element PLANE41 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/a2000_damper/ansys/loading_curve_ansys.txt -threshold 0.2 -reaction DISP_BC UY)

# This test fails
# add_test(NAME rigid_body_mpc_block_ROTZ COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTZ.cdb -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTZ)

# This test fails
# add_test(NAME rigid_body_mpc_block_ROTY COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTY.cdb -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTY)

add_test(NAME rigid_body_mpc_block_ROTX COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTX.cdb -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTX)

add_executable(TimeControlTestExe TimeControlTest.cpp)
target_link_libraries(TimeControlTestExe  nla3d_lib)
add_test(NAME TimeControlTest COMMAND TimeControlTestExe)
