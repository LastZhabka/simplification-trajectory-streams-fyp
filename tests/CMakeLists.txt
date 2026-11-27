add_executable(ssk_tests
  test_main.cpp
  test_trajectory.cpp
  test_simplify.cpp
  test_sweep.cpp
  test_frechet.cpp
  test_geometry.cpp
  test_json.cpp
  test_match.cpp
  test_measure.cpp
  test_greedy_stabbing.cpp
  test_regions.cpp
  test_bench.cpp
  test_corpus.cpp
  test_target.cpp
  test_size.cpp
)
target_link_libraries(ssk_tests PRIVATE ssk ssk_sweep ssk_measure ssk_bench_lib
                      ssk_target_lib ssk_size_lib)
target_include_directories(ssk_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME unit COMMAND ssk_tests)
