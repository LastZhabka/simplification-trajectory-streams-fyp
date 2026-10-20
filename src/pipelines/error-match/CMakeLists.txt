add_library(ssk_match_lib match.cpp)
target_link_libraries(ssk_match_lib PUBLIC ssk_measure)
target_include_directories(ssk_match_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ssk_match error_match.cpp)
target_link_libraries(ssk_match PRIVATE ssk_match_lib)

ssk_pipeline_warnings(ssk_match_lib ssk_match)
