cmake_minimum_required(VERSION 3.16)
project(wildfire CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Turn this OFF to build only the simulation core and the tests (no SFML or GDAL needed)
option(WF_BUILD_APP "Build the SFML/GDAL application" ON)

# Warnings for our own code only
function(wf_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra)
  endif()
endfunction()

# --- Simulation + planning: no SFML, no GDAL. This is what the tests link against. ---
add_library(wf_core STATIC
  src/distance_transform.cpp
  src/safe_planner.cpp
  src/voronoi.cpp
  src/boustrophedon.cpp
  src/aco.cpp
  src/mission_planner.cpp
  src/drone.cpp
  src/mission_controller.cpp
  src/fire_model.cpp
  src/fuel_models.cpp
  src/perimeter.cpp)
target_include_directories(wf_core PUBLIC src)
wf_warnings(wf_core)

# --- Headless tests ---
enable_testing()
add_executable(wf_tests tests/tests.cpp)
target_link_libraries(wf_tests PRIVATE wf_core)
wf_warnings(wf_tests)
add_test(NAME wf_tests COMMAND wf_tests)

if(WF_BUILD_APP)
  # The code uses the SFML 2 API (sf::Uint8, VideoMode(w, h), pollEvent(event)...), so SFML 3 won't do.
  find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)
  if(SFML_VERSION VERSION_GREATER_EQUAL 3)
    message(FATAL_ERROR "Found SFML ${SFML_VERSION}; this project needs SFML 2.5 - 2.6")
  endif()
  find_package(GDAL REQUIRED)

  # GIS input: the only code that touches GDAL
  add_library(wf_io STATIC src/raster_io.cpp src/terrain.cpp)
  target_link_libraries(wf_io PUBLIC wf_core GDAL::GDAL)
  wf_warnings(wf_io)

  # Drawing: the only code (besides main) that touches SFML
  add_library(wf_render STATIC src/renderer.cpp)
  target_link_libraries(wf_render PUBLIC wf_core sfml-graphics sfml-window sfml-system)
  wf_warnings(wf_render)

  add_executable(wildfire src/main.cpp)
  target_link_libraries(wildfire PRIVATE wf_core wf_io wf_render)
  wf_warnings(wildfire)

  # Run from the repo root (where elevation_map.tif lives) when launched from Visual Studio
  set_target_properties(wildfire PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()
