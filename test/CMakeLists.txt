# test divided by labels:
# FUNC - functional tests (unit tests) that checks some functionality works correct
# BENCH - benchmark tests that check how long and correct some FEM problem are solver


# FUNC tests
add_executable(mat_math_test_exe mat_math.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
target_link_libraries(mat_math_test_exe math)
add_test(NAME mat_math_test COMMAND mat_math_test_exe)
set_tests_properties(mat_math_test PROPERTIES LABELS "FUNC")


# TODO: undone! (see commented lines in the source code)
# add_executable(mat_sparse_pardiso_test_exe mat_sparse_pardiso.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
# target_link_libraries(mat_sparse_pardiso_test_exe math ${MKL_LIBRARIES})
# add_test(NAME mat_sparse_pardiso_test COMMAND mat_sparse_pardiso_test_exe)
# set_tests_properties(mat_sparse_pardiso_test PROPERTIES LABELS "FUNC")


add_executable(mat_prod_test_exe mat_prod/mat_prod_test.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
target_link_libraries(mat_prod_test_exe math)
add_test(NAME mat_prod_test COMMAND mat_prod_test_exe -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)
set_tests_properties(mat_prod_test PROPERTIES LABELS "FUNC")


add_executable(eigen_mat_prod_test_exe mat_prod/eigen_mat_prod_test.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
target_link_libraries(eigen_mat_prod_test_exe math)
add_test(NAME eigen_mat_prod_test COMMAND eigen_mat_prod_test_exe -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)
set_tests_properties(eigen_mat_prod_test PROPERTIES LABELS "FUNC")


# This test fails
# add_test(NAME rigid_body_mpc_block_ROTZ COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTZ.cdb -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTZ)
# set_tests_properties(rigid_body_mpc_block_ROTZ PROPERTIES LABELS "FUNC")


# This test fails
# add_test(NAME rigid_body_mpc_block_ROTY COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTY.cdb -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTY)
# set_tests_properties(rigid_body_mpc_block_ROTY PROPERTIES LABELS "FUNC")


add_test(NAME rigid_body_mpc_block_ROTX COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTX.cdb
    -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk
    -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt
    -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTX)
set_tests_properties(rigid_body_mpc_block_ROTX PROPERTIES LABELS "FUNC")


add_executable(TimeControlTestExe TimeControlTest.cpp)
target_link_libraries(TimeControlTestExe  nla3d_lib)
add_test(NAME TimeControlTest COMMAND TimeControlTestExe)
set_tests_properties(TimeControlTest PROPERTIES LABELS "FUNC")


add_test(NAME truss3D_test COMMAND nla3d_truss)
set_tests_properties(truss3D_test PROPERTIES LABELS "FUNC")


# BENCH tests
# this test takes ~90 sec on Release..
add_test(NAME 3d_damper COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/3d_damper/model.cdb
    -material Neo-Hookean 10 5000 -loadsteps 20 -novtk
    -refcurve ${PROJECT_SOURCE_DIR}/test/3d_damper/ansys/loading_curve_ansys.txt
    -reaction TOP_SIDE UY -threshold 3.0 -element SOLID81)
 set_tests_properties(3d_damper PROPERTIES LABELS "BENCH")
 

# this test takes toooo long with gauss eq. solver
add_test(NAME a2000_damper COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/a2000_damper/a2000.cdb
    -element PLANE41 -material Neo-Hookean 1 500 -loadsteps 20 -novtk
    -refcurve ${PROJECT_SOURCE_DIR}/test/a2000_damper/ansys/loading_curve_ansys.txt
    -threshold 0.2 -reaction DISP_BC UY)
set_tests_properties(a2000_damper PROPERTIES LABELS "BENCH")
