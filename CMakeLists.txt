cmake_minimum_required(VERSION 3.22)

find_program(CCACHE_PROGRAM ccache)

# If ccache is found, configure CMake to use it as a compiler launcher
if(CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    # Add other languages if needed, e.g., CMAKE_CUDA_COMPILER_LAUNCHER
    message(STATUS "Ccache found and enabled: ${CCACHE_PROGRAM}")
else()
    message(STATUS "Ccache not found. Compilation will proceed without caching.")
endif()


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
option(ENABLE_CUDA "Enable CUDA support" OFF)
if (ENABLE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)
endif()
# Set compiler options
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
   if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      # add_compile_options(-g -Wall -Wextra -O0 -Wpedantic -Wshadow -Wconversion -fsanitize=address,undefined )
      add_compile_options(-g -Wall -Wextra -O0)
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
get_target_property(cns_sources AMReX::amrex SOURCES)

if (ENABLE_CUDA)
    # set_source_files_properties(${cns_sources} PROPERTIES LANGUAGE CUDA)
    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        message(STATUS "No CUDA architectures specified. Native will be used.")
        set(CMAKE_CUDA_ARCHITECTURES native)
    endif()
    set_target_properties(${EXEC_NAME} PROPERTIES
                          CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES}
                          CUDA_SEPARABLE_COMPILATION ON)
    target_compile_features(${EXEC_NAME} PRIVATE cuda_std_17)
    target_compile_options(${EXEC_NAME} PRIVATE
                           $<$<COMPILE_LANGUAGE:CUDA>:-m64>
                           $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
                           $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>
                           $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>)
    target_link_libraries(${EXEC_NAME} PRIVATE CUDA::cudart)
endif()


if(ENABLE_TESTING)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()

install(TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_BINARY_DIR})
install(TARGETS "cns" DESTINATION ${CMAKE_BINARY_DIR})



