add_executable(ssk_tests
  test_main.cpp
  test_trajectory.cpp
)
target_link_libraries(ssk_tests PRIVATE ssk)
target_include_directories(ssk_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME unit COMMAND ssk_tests)
