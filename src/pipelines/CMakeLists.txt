# Each pipeline is a directory: the logic in a library so the tests can drive it, and a thin
# main beside it.
function(ssk_pipeline_warnings)
  foreach(target ${ARGV})
    if(MSVC)
      target_compile_options(${target} PRIVATE /W4)
    else()
      target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion)
    endif()
    if(SSK_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE $<IF:$<BOOL:${MSVC}>,/WX,-Werror>)
    endif()
  endforeach()
endfunction()

add_subdirectory(simplify-sweep)
add_subdirectory(frechet-distance)
