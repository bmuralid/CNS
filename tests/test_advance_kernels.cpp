#include <iostream>
#include <cmath>
#include <array>
#include <AMReX_Array4.H>
#include <AMReX_FArrayBox.H>
#include <AMReX.H>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_session.hpp>
#include "pyro.H"
#include "Thermo.H"
#include "Kernels.H"
#include "CNS_index_macros.H"

// Test tolerance for floating point comparisons
const double TOLERANCE = 1e-10;
using namespace amrex;

// AMReX setup and teardown fixtures
struct AMReXFixture {
    AMReXFixture() {
        int argc=0;
        char** argv = {nullptr};
        amrex::Initialize(argc, argv);
    }

    ~AMReXFixture() {
        amrex::Finalize();
    }
};

// Test function for derivative kernel - basic functionality
TEST_CASE_METHOD(AMReXFixture, "derivative_basic", "[advance]") {
    std::cout << "Testing derivative kernel basic functionality..." << std::endl;

    // Create test data - use a larger box to have neighbors
    const int i = 1, j = 1, k = 1;  // Center point
    const int ncomp_prims = NPRIM;
    const auto dx = amrex::Real(0.01);  // 1 cm grid spacing

    // Create FArrayBox for primitive variables (3x3x3 box for neighbors)
    amrex::Box prims_box(amrex::IntVect(0,0,0), amrex::IntVect(2,2,2));
    amrex::FArrayBox prims_fab(prims_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Create FArrayBox for derivatives
    amrex::FArrayBox deriv_fab(prims_box, ncomp_prims);
    auto deriv = deriv_fab.array();

    // Initialize primitive variables with a known function
    // Use a simple linear function: f(x) = 2*x + 1 for density
    // This should give a derivative of 2.0
    for (int ii = 0; ii <= 2; ++ii) {
        for (int jj = 0; jj <= 2; ++jj) {
            for (int kk = 0; kk <= 2; ++kk) {
                const auto x = amrex::Real(ii) * dx;  // x-coordinate
                prims(ii, jj, kk, QRHO) = 2.0 * x + 1.0;  // Linear in x
                prims(ii, jj, kk, QU) = 0.0;  // Zero velocity
                prims(ii, jj, kk, QV) = 0.0;
#if (AMREX_SPACEDIM == 3)
                prims(ii, jj, kk, QW) = 0.0;
#endif
                prims(ii, jj, kk, QPRES) = 101325.0;  // Constant pressure
                prims(ii, jj, kk, QTEMP) = 300.0;     // Constant temperature
                prims(ii, jj, kk, QY1) = 0.767;       // Constant mass fractions
                prims(ii, jj, kk, QY2) = 0.233;
            }
        }
    }

    // Test derivative in x-direction (idir = 0)
    int idir = 0;
    derivative(i, j, k, idir, 1, prims, deriv, dx);

    // Verify results
    // For f(x) = 2*x + 1, the derivative should be 2.0
    // But we need to account for the grid spacing: df/dx = 2.0, so df = 2.0 * dx
    amrex::Real expected_deriv = 2.0;  // This is df/dx, not df
    REQUIRE(std::abs(deriv(i, j, k, QRHO) - expected_deriv) < TOLERANCE);

    // Other components should be zero (constant functions)
    REQUIRE(std::abs(deriv(i, j, k, QU)) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QV)) < TOLERANCE);
#if (AMREX_SPACEDIM == 3)
    REQUIRE(std::abs(deriv(i, j, k, QW)) < TOLERANCE);
#endif
    REQUIRE(std::abs(deriv(i, j, k, QPRES)) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QTEMP)) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QY1)) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QY2)) < TOLERANCE);

    std::cout << "  ✓ derivative basic test passed" << std::endl;
    std::cout << "    Density derivative: " << deriv(i, j, k, QRHO) << " (expected: " << expected_deriv << ")" << std::endl;
}

// Test function for derivative kernel - all directions
TEST_CASE_METHOD(AMReXFixture, "derivative_directions", "[advance]") {
    std::cout << "Testing derivative kernel in all directions..." << std::endl;

    const int i = 1, j = 1, k = 1;
    const int ncomp_prims = NPRIM;
    const amrex::Real dx = 0.01;

    // Create FArrayBox for primitive variables
    amrex::Box prims_box(amrex::IntVect(0,0,0), amrex::IntVect(2,2,2));
    amrex::FArrayBox prims_fab(prims_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Create FArrayBox for derivatives
    amrex::FArrayBox deriv_fab(prims_box, ncomp_prims);
    auto deriv = deriv_fab.array();

    // Initialize with functions that vary in different directions
    for (int ii = 0; ii <= 2; ++ii) {
        for (int jj = 0; jj <= 2; ++jj) {
            for (int kk = 0; kk <= 2; ++kk) {
                const auto x = amrex::Real(ii) * dx;
                const auto y = amrex::Real(jj) * dx;
                const auto z = amrex::Real(kk) * dx;

                prims(ii, jj, kk, QRHO) = 1.0 + 0.5 * x;           // Linear in x
                prims(ii, jj, kk, QU) = 0.0;
                prims(ii, jj, kk, QV) = 0.0;
#if (AMREX_SPACEDIM == 3)
                prims(ii, jj, kk, QW) = 0.0;
#endif
                prims(ii, jj, kk, QPRES) = 101325.0 + 1000.0 * y;  // Linear in y
                prims(ii, jj, kk, QTEMP) = 300.0;
#if (AMREX_SPACEDIM == 3)
                prims(ii, jj, kk, QTEMP) = 300.0 + 50.0 * z;       // Linear in z (3D only)
#endif
                prims(ii, jj, kk, QY1) = 0.767;
                prims(ii, jj, kk, QY2) = 0.233;
            }
        }
    }

    // Test x-direction derivative
    derivative(i, j, k, 0, 1, prims, deriv, dx);
    REQUIRE(std::abs(deriv(i, j, k, QRHO) - 0.5) < TOLERANCE);  // df/dx = 0.5
    REQUIRE(std::abs(deriv(i, j, k, QPRES)) < TOLERANCE);        // Constant in x

    // Test y-direction derivative
    derivative(i, j, k, 1, 1, prims, deriv, dx);
    REQUIRE(std::abs(deriv(i, j, k, QPRES) - 1000.0) < TOLERANCE);  // df/dy = 1000.0
    REQUIRE(std::abs(deriv(i, j, k, QRHO)) < TOLERANCE);            // Constant in y

#if (AMREX_SPACEDIM == 3)
    // Test z-direction derivative (3D only)
    derivative(i, j, k, 2, 1, prims, deriv, dx);
    REQUIRE(std::abs(deriv(i, j, k, QTEMP) - 50.0) < TOLERANCE);    // df/dz = 50.0
    REQUIRE(std::abs(deriv(i, j, k, QRHO)) < TOLERANCE);            // Constant in z
    REQUIRE(std::abs(deriv(i, j, k, QPRES)) < TOLERANCE);           // Constant in z
#endif

    std::cout << "  ✓ derivative directions test passed" << std::endl;
}

// Test function for derivative kernel - edge cases
TEST_CASE_METHOD(AMReXFixture, "derivative_edge_cases", "[advance]") {
    std::cout << "Testing derivative kernel edge cases..." << std::endl;

    const int i = 1, j = 1, k = 1;
    const int ncomp_prims = NPRIM;
    const amrex::Real dx = 0.001;

    // Create FArrayBox for primitive variables
    amrex::Box prims_box(amrex::IntVect(0,0,0), amrex::IntVect(2,2,2));
    amrex::FArrayBox prims_fab(prims_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Create FArrayBox for derivatives
    amrex::FArrayBox deriv_fab(prims_box, ncomp_prims);
    auto deriv = deriv_fab.array();

    // Test case 1: Constant function (should give zero derivative)
    for (int ii = 0; ii <= 2; ++ii) {
        for (int jj = 0; jj <= 2; ++jj) {
            for (int kk = 0; kk <= 2; ++kk) {
                prims(ii, jj, kk, QRHO) = 1.2;  // Constant
                prims(ii, jj, kk, QU) = 0.0;
                prims(ii, jj, kk, QV) = 0.0;
#if (AMREX_SPACEDIM == 3)
                prims(ii, jj, kk, QW) = 0.0;
#endif
                prims(ii, jj, kk, QPRES) = 101325.0;
                prims(ii, jj, kk, QTEMP) = 300.0;
                prims(ii, jj, kk, QY1) = 0.767;
                prims(ii, jj, kk, QY2) = 0.233;
            }
        }
    }

    derivative(i, j, k, 0, 1, prims, deriv, dx);
    REQUIRE(std::abs(deriv(i, j, k, QRHO)) < TOLERANCE);  // Should be zero

    // Test case 2: Quadratic function (should give linear derivative)
    // f(x) = x^2, so df/dx = 2x
    for (int ii = 0; ii <= 2; ++ii) {
        for (int jj = 0; jj <= 2; ++jj) {
            for (int kk = 0; kk <= 2; ++kk) {
                const auto x = amrex::Real(ii) * dx;
                prims(ii, jj, kk, QRHO) = x * x;  // Quadratic in x
                prims(ii, jj, kk, QU) = 0.0;
                prims(ii, jj, kk, QV) = 0.0;
#if (AMREX_SPACEDIM == 3)
                prims(ii, jj, kk, QW) = 0.0;
#endif
                prims(ii, jj, kk, QPRES) = 101325.0;
                prims(ii, jj, kk, QTEMP) = 300.0;
                prims(ii, jj, kk, QY1) = 0.767;
                prims(ii, jj, kk, QY2) = 0.233;
            }
        }
    }

    derivative(i, j, k, 0, 0, prims, deriv, dx);
    // For f(x) = x^2, at x=1, df/dx = 2*1 = 2
    amrex::Real expected_deriv_quad = 2.0*dx;
    std::cout << "    Quadratic function derivative computed: " << deriv(i, j, k, QRHO) << std::endl;
    REQUIRE(std::abs(deriv(i, j, k, QRHO) - expected_deriv_quad) < 1e-6);  // Allow some tolerance for finite difference

    std::cout << "  ✓ derivative edge cases test passed" << std::endl;
    std::cout << "    Constant function derivative: " << deriv(i, j, k, QRHO) << " (should be ~0)" << std::endl;
}

// Test function for derivative kernel - monocent limiter behavior
TEST_CASE_METHOD(AMReXFixture, "derivative_monocent", "[advance]") {
    std::cout << "Testing derivative kernel monocent limiter behavior..." << std::endl;

    const int i = 1, j = 1, k = 1;
    const int ncomp_prims = NPRIM;
    const amrex::Real dx = 1.0;

    // Create FArrayBox for primitive variables
    amrex::Box prims_box(amrex::IntVect(0,0,0), amrex::IntVect(2,2,2));
    amrex::FArrayBox prims_fab(prims_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Create FArrayBox for derivatives
    amrex::FArrayBox deriv_fab(prims_box, ncomp_prims);
    auto deriv = deriv_fab.array();

    // Test case: Function with different slopes on left and right
    // This tests the monocent limiter which should choose the most conservative estimate
    for (int ii = 0; ii <= 2; ++ii) {
        for (int jj = 0; jj <= 2; ++jj) {
            for (int kk = 0; kk <= 2; ++kk) {
                if (ii == 0) {
                    prims(ii, jj, kk, QRHO) = 0.0;  // Left value
                } else if (ii == 1) {
                    prims(ii, jj, kk, QRHO) = 1.0;  // Center value
                } else {
                    prims(ii, jj, kk, QRHO) = 3.0;  // Right value (steeper slope)
                }
                prims(ii, jj, kk, QU) = 0.0;
                prims(ii, jj, kk, QV) = 0.0;
#if (AMREX_SPACEDIM == 3)
                prims(ii, jj, kk, QW) = 0.0;
#endif
                prims(ii, jj, kk, QPRES) = 101325.0;
                prims(ii, jj, kk, QTEMP) = 300.0;
                prims(ii, jj, kk, QY1) = 0.767;
                prims(ii, jj, kk, QY2) = 0.233;
            }
        }
    }

    derivative(i, j, k, 0, 1, prims, deriv, dx);

    // The monocent limiter should choose the most conservative (smallest magnitude) derivative
    // Left slope: (1-0)/dx = 1/dx
    // Central slope: (3-0)/(2*dx) = 1.5/dx
    // Right slope: (3-1)/dx = 2/dx
    // Monocent should choose the minimum positive value: 1/dx
    amrex::Real expected_deriv_monocent = 1.0;  // This is df/dx
    REQUIRE(std::abs(deriv(i, j, k, QRHO) - expected_deriv_monocent) < TOLERANCE);

    std::cout << "  ✓ derivative monocent test passed" << std::endl;
    std::cout << "    Monocent limited derivative: " << deriv(i, j, k, QRHO) << " (expected: " << expected_deriv_monocent << ")" << std::endl;
}

// Test function for derivative kernel - grid spacing dependence
TEST_CASE_METHOD(AMReXFixture, "derivative_grid_spacing", "[advance]") {
    std::cout << "Testing derivative kernel grid spacing dependence..." << std::endl;

    const int i = 1, j = 1, k = 1;
    const int ncomp_prims = NPRIM;

    // Create FArrayBox for primitive variables
    amrex::Box prims_box(amrex::IntVect(0,0,0), amrex::IntVect(2,2,2));
    amrex::FArrayBox prims_fab(prims_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Create FArrayBox for derivatives
    amrex::FArrayBox deriv_fab(prims_box, ncomp_prims);
    auto deriv = deriv_fab.array();
    // Test with different grid spacings
    std::array<amrex::Real, 3> dx_values = {0.01, 0.005, 0.02};
    amrex::Real expected_deriv = 3.0;  // df/dx = 3.0

    for (auto dx : dx_values) {
      // Initialize with linear function f(x) = 3*x
      for (int ii = 0; ii <= 2; ++ii) {
          for (int jj = 0; jj <= 2; ++jj) {
              for (int kk = 0; kk <= 2; ++kk) {
                  const auto x = amrex::Real(ii) * dx;
                  prims(ii, jj, kk, QRHO) = 3.0 * x;  // Linear function
                  prims(ii, jj, kk, QU) = 0.0;
                  prims(ii, jj, kk, QV) = 0.0;
#if (AMREX_SPACEDIM == 3)
                  prims(ii, jj, kk, QW) = 0.0;
#endif
                  prims(ii, jj, kk, QPRES) = 101325.0;
                  prims(ii, jj, kk, QTEMP) = 300.0;
                  prims(ii, jj, kk, QY1) = 0.767;
                  prims(ii, jj, kk, QY2) = 0.233;
              }
          }
      }
      derivative(i, j, k, 0, 1, prims, deriv, dx);
      REQUIRE(std::abs(deriv(i, j, k, QRHO) - expected_deriv) < TOLERANCE);
    }

    std::cout << "  ✓ derivative grid spacing test passed" << std::endl;
    std::cout << "    Tested with dx = 0.01, 0.005, 0.02 - all gave derivative = " << deriv(i, j, k, QRHO) << std::endl;
}

// Test function for derivative kernel - all primitive variables
TEST_CASE_METHOD(AMReXFixture, "derivative_all_variables", "[advance]") {
    std::cout << "Testing derivative kernel for all primitive variables..." << std::endl;

    const int i = 1, j = 1, k = 1;
    const int ncomp_prims = NPRIM;
    const amrex::Real dx = 1.0;

    // Create FArrayBox for primitive variables
    amrex::Box prims_box(amrex::IntVect(0,0,0), amrex::IntVect(2,2,2));
    amrex::FArrayBox prims_fab(prims_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Create FArrayBox for derivatives
    amrex::FArrayBox deriv_fab(prims_box, ncomp_prims);
    auto deriv = deriv_fab.array();

    // Initialize all primitive variables with different linear functions
    for (int ii = 0; ii <= 2; ++ii) {
        for (int jj = 0; jj <= 2; ++jj) {
            for (int kk = 0; kk <= 2; ++kk) {
                const auto x = amrex::Real(ii) * dx;
                prims(ii, jj, kk, QRHO) = 1.0 + 0.5 * x;           // Linear in x
                prims(ii, jj, kk, QU) = 10.0 + 2.0 * x;           // Linear in x
                prims(ii, jj, kk, QV) = 5.0 + 1.0 * x;            // Linear in x
#if (AMREX_SPACEDIM == 3)
                prims(ii, jj, kk, QW) = 3.0 + 0.5 * x;            // Linear in x
#endif
                prims(ii, jj, kk, QPRES) = 101325.0 + 1000.0 * x; // Linear in x
                prims(ii, jj, kk, QTEMP) = 300.0 + 10.0 * x;     // Linear in x
                prims(ii, jj, kk, QY1) = 0.767 + 0.01 * x;       // Linear in x
                prims(ii, jj, kk, QY2) = 0.233 - 0.01 * x;       // Linear in x
            }
        }
    }

    derivative(i, j, k, 0, 1, prims, deriv, dx);

    // Verify all derivatives
    REQUIRE(std::abs(deriv(i, j, k, QRHO) - 0.5) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QU) - 2.0) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QV) - 1.0) < TOLERANCE);
#if (AMREX_SPACEDIM == 3)
    REQUIRE(std::abs(deriv(i, j, k, QW) - 0.5) < TOLERANCE);
#endif
    REQUIRE(std::abs(deriv(i, j, k, QPRES) - 1000.0) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QTEMP) - 10.0) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QY1) - 0.01) < TOLERANCE);
    REQUIRE(std::abs(deriv(i, j, k, QY2) - (-0.01)) < TOLERANCE);

    std::cout << "  ✓ derivative all variables test passed" << std::endl;
    std::cout << "    All " << NPRIM << " primitive variables tested successfully" << std::endl;
}

// Global setup and teardown for AMReX
struct AMReXGlobalSetup {
    AMReXGlobalSetup() {
        // std::cout << "Global AMReX setup..." << std::endl;
        // std::cout << "AMREX_SPACEDIM = " << AMREX_SPACEDIM << std::endl;
        // std::cout << "NCONS = " << NCONS << std::endl;
        // std::cout << "NPRIM = " << NPRIM << std::endl;
        // std::cout << "NSP = " << NSP << std::endl;
        // std::cout << std::endl;
        /* int argc=0; */
        /* char** argv = {nullptr}; */
        /* amrex::Initialize(argc, argv); */

    }

    ~AMReXGlobalSetup() {
        // std::cout << "Global AMReX teardown..." << std::endl;
        /* amrex::Finalize(); */
    }
};

// Global instance for setup/teardown
AMReXGlobalSetup global_setup;
