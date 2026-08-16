add_library(ssk_audit_lib audit.cpp)
target_link_libraries(ssk_audit_lib PUBLIC ssk)
target_include_directories(ssk_audit_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ssk_audit audit_corpus.cpp)
target_link_libraries(ssk_audit PRIVATE ssk_audit_lib)

ssk_pipeline_warnings(ssk_audit_lib ssk_audit)
