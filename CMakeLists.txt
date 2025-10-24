cmake_minimum_required(VERSION 3.22)

project(
    "amrex_cns"
    VERSION 0.0.0
    LANGUAGES C CXX Fortran)

# Global CMake variables are set here
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_TESTING "Enable testing" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
# Find AMREx library
find_package(AMReX REQUIRED HINTS $ENV{AMREX_ROOT_DIR})

# Find Catch2 for testing
if(ENABLE_TESTING)
    find_package(Catch2 REQUIRED HINTS $ENV{CATCH2_ROOT_DIR})
endif()

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

message(${CMAKE_CXX_COMPILER_ID})

option(USE_MPI "Enable MPI" ON)

# Set compiler options
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
   if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      # add_compile_options(-g -Wall -Wextra -O0 -Wpedantic -Wshadow -Wconversion -fsanitize=address,undefined )
      add_compile_options(-g -Wall -Wextra -O0 )
   elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
       add_compile_options(-Wall -Wextra -O2 -g)
   else()
      add_compile_options(-Wall -Wextra -O2 -g)
   endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "NVHPC")
   add_compile_options( -O2)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
   add_compile_options(-O2 -DIFORT)
endif()

if(${USE_MPI})
   find_package(MPI REQUIRED COMPONENTS)
   # include_directories(${MPI_INCLUDE_PATH})
   set(EXEC_NAME "${PROJECT_NAME}_${CMAKE_CXX_COMPILER_ID}.mpi")
else()
    set(EXEC_NAME "${PROJECT_NAME}_${CMAKE_CXX_COMPILER_ID}.serial")
endif()

add_subdirectory(src)
add_subdirectory(exec)

if(ENABLE_TESTING)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()

install(TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_BINARY_DIR})
install(TARGETS "cns" DESTINATION ${CMAKE_BINARY_DIR})



