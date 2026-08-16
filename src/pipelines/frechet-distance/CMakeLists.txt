add_library(ssk_measure measure.cpp)
target_link_libraries(ssk_measure PUBLIC ssk)
target_include_directories(ssk_measure PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ssk_frechet frechet_distances.cpp)
target_link_libraries(ssk_frechet PRIVATE ssk_measure)

ssk_pipeline_warnings(ssk_measure ssk_frechet)
