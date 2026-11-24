add_library(ssk_bench_lib bench.cpp)
target_link_libraries(ssk_bench_lib PUBLIC ssk)
target_include_directories(ssk_bench_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ssk_bench benchmark.cpp)
target_link_libraries(ssk_bench PRIVATE ssk_bench_lib ssk_match_lib)

ssk_pipeline_warnings(ssk_bench_lib ssk_bench)
