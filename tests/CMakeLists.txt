# CMakeLists.txt for CNS tests
include(Catch)

# Find AMReX library
# find_package(AMReX REQUIRED HINTS $ENV{AMREX_ROOT_DIR})

# Create test executables
add_executable(test_thermo_kernels test_thermo_kernels.cpp)
add_executable(test_advance_kernels test_advance_kernels.cpp)

# Link with AMReX and Catch2
target_link_libraries(test_thermo_kernels PUBLIC AMReX::amrex ${EXEC_NAME})
target_link_libraries(test_advance_kernels PUBLIC AMReX::amrex ${EXEC_NAME})
# target_link_libraries(test_thermo_kernels PUBLIC ${EXEC_NAME})
# target_link_libraries(test_advance_kernels PUBLIC ${EXEC_NAME})
target_link_libraries(test_thermo_kernels PRIVATE Catch2::Catch2WithMain)
target_link_libraries(test_advance_kernels PRIVATE Catch2::Catch2WithMain)


# Include directories
target_include_directories(test_thermo_kernels PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${AMREX_INCLUDE_DIRS}
)

target_include_directories(test_advance_kernels PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${AMREX_INCLUDE_DIRS}
)

# Add compile definitions
target_compile_definitions(test_thermo_kernels PRIVATE
    AMREX_SPACEDIM=${AMREX_SPACEDIM}
)

target_compile_definitions(test_advance_kernels PRIVATE
    AMREX_SPACEDIM=${AMREX_SPACEDIM}
)

# Set C++ standard
set_target_properties(test_thermo_kernels PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

set_target_properties(test_advance_kernels PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

# Add test to CTest
# add_test(NAME thermo_kernels_test COMMAND test_thermo_kernels)
# add_test(NAME advance_kernels_test COMMAND test_advance_kernels)
catch_discover_tests("test_thermo_kernels")
catch_discover_tests("test_advance_kernels")
