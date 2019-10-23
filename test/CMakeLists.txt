add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${PROJECT_SOURCE_DIR}/test)

add_library(test_main main.cpp)
target_link_libraries(test_main Catch)

set(sources
  # odes
  test_harmonic_oscillator_ode.cpp)


foreach (sourcefile ${sources})
  string(REPLACE ".cpp" "" sourcename ${sourcefile})
  add_executable(${sourcename} ${sourcefile})
  #target_link_libraries(${sourcename} sam)
  target_link_libraries(${sourcename} test_main)
  add_test(${sourcename} ${sourcename})
endforeach(sourcefile ${sources})
