add_library(ssk_sweep sweep.cpp)
target_link_libraries(ssk_sweep PUBLIC ssk)
target_include_directories(ssk_sweep PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ssk_simplify simplify_sweep.cpp)
target_link_libraries(ssk_simplify PRIVATE ssk_sweep)

ssk_pipeline_warnings(ssk_sweep ssk_simplify)
