# CMakeLists.txt for CNS tests

# Find AMReX library
find_package(AMReX REQUIRED HINTS $ENV{AMREX_ROOT_DIR})

# Create test executable
add_executable(test_thermo_kernels test_thermo_kernels.cpp)

# Link with AMReX
target_link_libraries(test_thermo_kernels AMReX::amrex ${EXEC_NAME})


# Include directories
target_include_directories(test_thermo_kernels PRIVATE 
    ${CMAKE_SOURCE_DIR}/src
    ${AMREX_INCLUDE_DIRS}
)

# Add compile definitions
target_compile_definitions(test_thermo_kernels PRIVATE
    AMREX_SPACEDIM=${AMREX_SPACEDIM}
)

# Set C++ standard
set_target_properties(test_thermo_kernels PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

# Add test to CTest
enable_testing()
add_test(NAME thermo_kernels_test COMMAND test_thermo_kernels)
