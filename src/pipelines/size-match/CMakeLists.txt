add_library(ssk_size_lib size.cpp)
target_link_libraries(ssk_size_lib PUBLIC ssk_measure ssk_match_lib)
target_include_directories(ssk_size_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ssk_size_match size_match.cpp)
target_link_libraries(ssk_size_match PRIVATE ssk_size_lib ssk_corpus)

ssk_pipeline_warnings(ssk_size_lib ssk_size_match)
