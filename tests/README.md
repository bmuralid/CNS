# CNS Thermo Kernels Tests

This directory contains tests for the CNS thermodynamic kernel functions:

- `cns_ctoprim`: Converts conservative variables to primitive variables
- `cns_primtoc`: Converts primitive variables to conservative variables  
- `cns_computedt`: Computes the time step based on CFL condition

## Test Functions

### `test_cns_ctoprim()`
Tests the conversion from conservative to primitive variables:
- Sets up conservative variables (density, momentum, energy, species)
- Calls `cns_ctoprim` to convert to primitive variables
- Verifies that density, velocity components, pressure, and temperature are correctly computed
- Checks that all values are physically reasonable (positive pressure, temperature, etc.)

### `test_cns_primtoc()`
Tests the conversion from primitive to conservative variables:
- Sets up primitive variables (density, velocity, pressure, temperature, species)
- Calls `cns_primtoc` to convert to conservative variables
- Verifies that density, momentum, total energy, and species mass are correctly computed
- Checks that all values are physically reasonable

### `test_cns_computedt()`
Tests the time step calculation:
- Sets up primitive variables with known velocity and thermodynamic state
- Calls `cns_computedt` with a given grid spacing
- Verifies that the computed time step is positive and reasonable
- Checks that the time step follows CFL condition: dt = dx / (|u| + c)

### `test_roundtrip()`
Tests that the conversions are inverse operations:
- Starts with conservative variables
- Converts to primitive variables using `cns_ctoprim`
- Converts back to conservative variables using `cns_primtoc`
- Verifies that the final conservative variables match the original (within numerical tolerance)

### `test_edge_cases()`
Tests edge cases:
- Very small density values
- Ensures functions handle edge cases gracefully without crashing

## Building and Running

### Using CMake (Recommended)
```bash
cd /path/to/CNS
mkdir build && cd build
cmake ..
make test_thermo_kernels
./test_thermo_kernels
```

### Using Makefile
```bash
cd tests
make
make test
```

### Using CTest
```bash
cd build
ctest -R thermo_kernels_test
```

## Test Data

The tests use realistic fluid dynamics data:
- **Density**: 1.2 kg/m³ (typical for air at standard conditions)
- **Velocity**: 100 m/s in x-direction, 50 m/s in y-direction, 25 m/s in z-direction (3D)
- **Pressure**: 101325 Pa (1 atmosphere)
- **Temperature**: 300 K (room temperature)
- **Species**: 76.7% O₂, 23.3% N₂ (air composition)
- **Grid spacing**: 0.01 m (1 cm)

## Expected Output

The test should output:
```
Starting CNS thermo kernels tests...
AMREX_SPACEDIM = 3
NCONS = 7
NPRIM = 8
NSP = 2

Testing cns_ctoprim...
  ✓ cns_ctoprim test passed
    Density: 1.2
    Velocity x: 100
    Velocity y: 50
    Velocity z: 25
    Pressure: [computed value]
    Temperature: [computed value]

Testing cns_primtoc...
  ✓ cns_primtoc test passed
    [similar output]

Testing cns_computedt...
  ✓ cns_computedt test passed
    Computed dt: [computed value]
    Expected order of magnitude: [expected value]

Testing round-trip conversion (cons -> prims -> cons)...
  ✓ Round-trip test passed
    [comparison values]

Testing edge cases...
  ✓ Edge case test passed

All tests passed successfully!
```

## Dependencies

- AMReX library
- C++20 compatible compiler
- The CNS source code (pyro.H, Thermo.H, Kernels.H, CNS_index_macros.H)

## Notes

- The tests use `assert()` statements for validation
- Floating-point comparisons use appropriate tolerances
- The tests are designed to be run in both 2D and 3D configurations
- All tests should pass for the functions to be considered working correctly
