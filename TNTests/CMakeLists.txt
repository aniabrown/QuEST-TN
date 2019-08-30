# Define all TN tests. These are currently all integration tests

set(test_names
    testCircuit1
    testCircuit1BLAS
    testCircuit2
    testCircuit3
    )

foreach(test_name ${test_names})
    # For all test circuits, add a building and running test
    add_executable(bin_${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.c)
    target_link_libraries(bin_${test_name} TN m)
    add_test(build_${test_name} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target bin_${test_name})
    add_test(test_${test_name} bin_${test_name})
endforeach()




