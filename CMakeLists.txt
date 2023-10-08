CMAKE_MINIMUM_REQUIRED(VERSION 3.13)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(simulation_cpp)

# Submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/eigen/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

##########################################
####### CPP EXECUTABLE COMPILATION #######
##########################################

set(SOURCE_FILES
  src/spring.cpp
  src/clock.cpp
  src/gravity.cpp
  src/two_body_helper.cpp
  src/mass_spring.cpp
  src/edge.cpp
  src/mesh_boundary.cpp
  src/integrate.cpp
  src/utility_functions.cpp
  src/gen_mesh_plane.cpp
  src/contact.cpp
)

set(SOURCE_FILES_MAIN
  ${SOURCE_FILES}
  src/main.cpp
  )

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/raylib/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/)


# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wall)
# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -g)
# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -O3)

add_executable(${PROJECT_NAME} ${SOURCE_FILES_MAIN})

target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/raylib/src/
)

target_link_libraries(${PROJECT_NAME} PUBLIC
  raylib
)

##########################################
####### PYTHON LIBRARY COMPILATION #######
##########################################

set(PYTHON_MODULE_NAME "symulathon")

set(PYTHON_MODULE_SRC
  ${SOURCE_FILES}
  src/pysimulation.cpp
  src/pybinding.cpp
  )

pybind11_add_module(${PYTHON_MODULE_NAME} MODULE ${PYTHON_MODULE_SRC})

target_include_directories(${PYTHON_MODULE_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/pybind/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/
  ${RENDERER_DIR}/include/
)

target_link_libraries(${PYTHON_MODULE_NAME} PUBLIC
  raylib
)
