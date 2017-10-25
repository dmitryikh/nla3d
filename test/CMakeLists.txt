# test divided by labels:
# FUNC - functional tests (unit tests) that checks some functionality works correct
# BENCH - benchmark tests that check how long and correct some FEM problem are solver


# FUNC tests
# ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp
set (TEST_SOURCES "mat_math.cpp")
set (TEST_NAME "mat_math")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})


# TODO: undone! (see commented lines in the source code)
# add_executable(mat_sparse_pardiso_test_exe mat_sparse_pardiso.cpp ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)
# target_link_libraries(mat_sparse_pardiso_test_exe math ${MKL_LIBRARIES})
# add_test(NAME mat_sparse_pardiso_test COMMAND mat_sparse_pardiso_test_exe)
# set_tests_properties(mat_sparse_pardiso_test PROPERTIES LABELS "FUNC")


set (TEST_SOURCES "mat_prod_test.cpp")
set (TEST_NAME "mat_prod")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})


set (TEST_SOURCES "eigen_mat_prod_test.cpp")
set (TEST_NAME "eigen_mat_prod")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} -dir ${PROJECT_SOURCE_DIR}/test/mat_prod/matlab -num 100)
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})


# This test fails
# add_test(NAME rigid_body_mpc_block_ROTZ COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTZ.cdb -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTZ)
# set_tests_properties(rigid_body_mpc_block_ROTZ PROPERTIES LABELS "FUNC")


# This test fails
# add_test(NAME rigid_body_mpc_block_ROTY COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTY.cdb -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTY)
# set_tests_properties(rigid_body_mpc_block_ROTY PROPERTIES LABELS "FUNC")


set (TEST_NAME "rigid_body_mpc_block_ROTX")
add_test(NAME ${TEST_NAME} COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/block_ROTX.cdb
    -element SOLID81 -material Neo-Hookean 1 500 -loadsteps 20 -novtk
    -refcurve ${PROJECT_SOURCE_DIR}/test/rigid_body_mpc/reference_MOMZ_reaction.txt
    -threshold 0.0001 -rigidbody 9 TOP_SIDE -reaction MASTER_NODE ROTX)
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check nla3d)


set (TEST_SOURCES "TimeControlTest.cpp")
set (TEST_NAME "TimeControl")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})


set (TEST_SOURCES "isoparametric_test.cpp")
set (TEST_NAME "Isoparametric")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})


set (TEST_NAME "truss3D")
add_test(NAME ${TEST_NAME} COMMAND nla3d_truss)
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check nla3d_truss)


set (TEST_SOURCES "sparse_math.cpp")
set (TEST_NAME "SparseMath")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})


set (TEST_SOURCES "QUADTH_test.cpp;${PROJECT_SOURCE_DIR}/src/FEReaders.cpp")
set (TEST_NAME "QUADTH_test")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${CMAKE_SOURCE_DIR}/test/QUADTH/heat.cdb
  ${CMAKE_SOURCE_DIR}/test/QUADTH/temp.txt)
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})


set (TEST_SOURCES "QUADTH_transient.cpp;${PROJECT_SOURCE_DIR}/src/FEReaders.cpp")
set (TEST_NAME "QUADTH_transient")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${CMAKE_SOURCE_DIR}/test/QUADTH/heat.cdb
  ${CMAKE_SOURCE_DIR}/test/QUADTH/temphistory.txt)
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})

set (TEST_SOURCES "fereaders.cpp;${PROJECT_SOURCE_DIR}/src/FEReaders.cpp")
set (TEST_NAME "fereaders")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${CMAKE_SOURCE_DIR}/test/QUADTH/heat.cdb)
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})

set (TEST_SOURCES "TETRA0_test.cpp;${PROJECT_SOURCE_DIR}/src/FEReaders.cpp")
set (TEST_NAME "TETRA0_test")
add_executable(${TEST_NAME} ${TEST_SOURCES})
target_link_libraries(${TEST_NAME} nla3d_lib)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${CMAKE_SOURCE_DIR}/test/TETRA0/tetra.cdb
    ${CMAKE_SOURCE_DIR}/test/TETRA0/disp.res ${CMAKE_SOURCE_DIR}/test/TETRA0/stress.res)
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "FUNC")
add_dependencies(check ${TEST_NAME})

# BENCH tests

# this test takes toooo long with gauss eq. solver
add_test(NAME a2000_damper COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/a2000_damper/a2000.cdb
    -element PLANE41 -material Neo-Hookean 1 500 -loadsteps 20 -novtk
    -refcurve ${PROJECT_SOURCE_DIR}/test/a2000_damper/ansys/loading_curve_ansys.txt
    -threshold 0.2 -reaction DISP_BC UY)
set_tests_properties(a2000_damper PROPERTIES LABELS "BENCH")

# this test takes ~90 sec on Release..
add_test(NAME 3d_damper COMMAND nla3d ${PROJECT_SOURCE_DIR}/test/3d_damper/model.cdb
  -material Neo-Hookean 10 5000 -loadsteps 20 -novtk
  -refcurve ${PROJECT_SOURCE_DIR}/test/3d_damper/ansys/loading_curve_ansys.txt
  -reaction TOP_SIDE UY -threshold 3.0 -element SOLID81)
set_tests_properties(3d_damper PROPERTIES LABELS "BENCH")
 
