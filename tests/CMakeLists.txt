add_executable(ssk_tests
  test_main.cpp
  test_trajectory.cpp
  test_simplify.cpp
  test_sweep.cpp
  test_frechet.cpp
  test_measure.cpp
  test_audit.cpp
)
target_link_libraries(ssk_tests PRIVATE ssk ssk_sweep ssk_measure ssk_audit_lib)
target_include_directories(ssk_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME unit COMMAND ssk_tests)
