add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${PROJECT_SOURCE_DIR}/test)

add_library(test_main main.cpp)
target_link_libraries(test_main Catch sam)

set(sources
  #systems
  test_generic_system.cpp
  test_rk4_system.cpp
  # odes
  test_harmonic_oscillator_ode.cpp
  #observers
  test_position_observer.cpp
  # options
  test_options.cpp)


foreach (sourcefile ${sources})
  string(REPLACE ".cpp" "" sourcename ${sourcefile})
  add_executable(${sourcename} ${sourcefile})
  #target_link_libraries(${sourcename} sam)
  target_link_libraries(${sourcename} test_main ${Boost_LIBRARIES})
  add_test(${sourcename} ${sourcename})
endforeach(sourcefile ${sources})
