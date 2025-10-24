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

// Test function for cns_ctoprim
TEST_CASE_METHOD(AMReXFixture, "cns_ctoprim", "[thermo]") {
    std::cout << "Testing cns_ctoprim..." << std::endl;

    /* int argc=0; */
    /* char** argv = {nullptr}; */
    /* amrex::Initialize(argc, argv); */
    // Create test data
    const int i = 0, j = 0, k = 0;
    const int ncomp_cons = NCONS;
    const int ncomp_prims = NPRIM;

    // Create FArrayBox for conservative variables
    amrex::Box cons_box(amrex::IntVect(0, 0, 0), amrex::IntVect(1, 1, 1));

    amrex::FArrayBox cons_fab(cons_box, ncomp_cons);
    auto cons = cons_fab.array();

    // Create FArrayBox for primitive variables
    amrex::FArrayBox prims_fab(cons_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Initialize conservative variables with test data
    // Density
    cons(i, j, k, URHO) = 1.2;  // kg/m^3

    // Momentum components
    cons(i, j, k, UMX) = 1.2 * 100.0;  // rho * ux
    cons(i, j, k, UMY) = 1.2 * 50.0;   // rho * uy
#if (AMREX_SPACEDIM == 3)
    cons(i, j, k, UMZ) = 1.2 * 25.0;   // rho * uz
#endif
    // Total energy (internal + kinetic)
    double ux = 100.0, uy = 50.0;
#if (AMREX_SPACEDIM == 3)
    double uz = 25.0;
#else
    double uz = 0.0;
#endif
    double kinetic_energy = 0.5 * 1.2 * (ux*ux + uy*uy + uz*uz);
    double internal_energy = 200000.0;  // J/kg
    cons(i, j, k, UEINT) = 1.2 * internal_energy + kinetic_energy;

    // Species mass fractions (stored as rho * Y)
    cons(i, j, k, UY1) = 1.2 * 0.767;  // O2 mass fraction
    cons(i, j, k, UY2) = 1.2 * 0.233;  // N2 mass fraction

    // Initialize primitive variables with initial guess
    prims(i, j, k, QTEMP) = 300.0;  // Initial temperature guess

    // Create thermo object
    pyro::pyro<double> thermo;

    // Call the function
    cns_ctoprim(i, j, k, cons, prims, thermo);

    // Verify results
    REQUIRE(std::abs(prims(i, j, k, QRHO) - 1.2) < TOLERANCE);
    REQUIRE(std::abs(prims(i, j, k, QU) - 100.0) < TOLERANCE);
    REQUIRE(std::abs(prims(i, j, k, QV) - 50.0) < TOLERANCE);
#if (AMREX_SPACEDIM == 3)
    REQUIRE(std::abs(prims(i, j, k, QW) - 25.0) < TOLERANCE);
#endif

    // Check that pressure and temperature are positive
    REQUIRE(prims(i, j, k, QPRES) > 0.0);
    REQUIRE(prims(i, j, k, QTEMP) > 0.0);

    std::cout << "  ✓ cns_ctoprim test passed" << std::endl;
    std::cout << "    Density: " << prims(i, j, k, QRHO) << std::endl;
    std::cout << "    Velocity x: " << prims(i, j, k, QU) << std::endl;
    std::cout << "    Velocity y: " << prims(i, j, k, QV) << std::endl;
#if (AMREX_SPACEDIM == 3)
    std::cout << "    Velocity z: " << prims(i, j, k, QW) << std::endl;
#endif
    std::cout << "    Pressure: " << prims(i, j, k, QPRES) << std::endl;
    std::cout << "    Temperature: " << prims(i, j, k, QTEMP) << std::endl;
    /* amrex::Finalize(); */
}

// Test function for cns_primtoc
TEST_CASE_METHOD(AMReXFixture, "cns_primtoc", "[thermo]") {
    std::cout << "Testing cns_primtoc..." << std::endl;

    // Create test data
    const int i = 0, j = 0, k = 0;
    const int ncomp_cons = NCONS;
    const int ncomp_prims = NPRIM;

    // Create FArrayBox for primitive variables
    amrex::Box prims_box(amrex::IntVect(0,0,0), amrex::IntVect(1,1,1));
    amrex::FArrayBox prims_fab(prims_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Create FArrayBox for conservative variables
    amrex::FArrayBox cons_fab(prims_box, ncomp_cons);
    auto cons = cons_fab.array();

    // Initialize primitive variables with test data
    prims(i, j, k, QRHO) = 1.2;      // kg/m^3
    prims(i, j, k, QU) = 100.0;      // m/s
    prims(i, j, k, QV) = 50.0;       // m/s
#if (AMREX_SPACEDIM == 3)
    prims(i, j, k, QW) = 25.0;       // m/s
#endif
    prims(i, j, k, QPRES) = 101325.0; // Pa (1 atm)
    prims(i, j, k, QTEMP) = 300.0;   // K

    // Species mass fractions
    prims(i, j, k, QY1) = 0.767;     // O2 mass fraction
    prims(i, j, k, QY2) = 0.233;     // N2 mass fraction

    // Create thermo object
    pyro::pyro<double> thermo;

    // Call the function
    cns_primtoc(i, j, k, prims, cons, thermo);

    // Verify results
    REQUIRE(std::abs(cons(i, j, k, URHO) - 1.2) < TOLERANCE);
    REQUIRE(std::abs(cons(i, j, k, UMX) - 1.2 * 100.0) < TOLERANCE);
    REQUIRE(std::abs(cons(i, j, k, UMY) - 1.2 * 50.0) < TOLERANCE);
#if (AMREX_SPACEDIM == 3)
    REQUIRE(std::abs(cons(i, j, k, UMZ) - 1.2 * 25.0) < TOLERANCE);
#endif

    // Check that total energy is positive
    /* REQUIRE(cons(i, j, k, UEINT) > 0.0); */

    // Check species mass fractions
    REQUIRE(std::abs(cons(i, j, k, UY1) - 1.2 * 0.767) < TOLERANCE);
    REQUIRE(std::abs(cons(i, j, k, UY2) - 1.2 * 0.233) < TOLERANCE);

    std::cout << "  ✓ cns_primtoc test passed" << std::endl;
    std::cout << "    Density: " << cons(i, j, k, URHO) << std::endl;
    std::cout << "    Momentum x: " << cons(i, j, k, UMX) << std::endl;
    std::cout << "    Momentum y: " << cons(i, j, k, UMY) << std::endl;
#if (AMREX_SPACEDIM == 3)
    std::cout << "    Momentum z: " << cons(i, j, k, UMZ) << std::endl;
#endif
    std::cout << "    Total energy: " << cons(i, j, k, UEINT) << std::endl;
    std::cout << "    Species 1 mass: " << cons(i, j, k, UY1) << std::endl;
    std::cout << "    Species 2 mass: " << cons(i, j, k, UY2) << std::endl;
}

// Test function for cns_computedt
TEST_CASE_METHOD(AMReXFixture, "cns_computedt", "[thermo]") {
    std::cout << "Testing cns_computedt..." << std::endl;

    // Create test data
    const int i = 0, j = 0, k = 0;
    const int ncomp_prims = NPRIM;

    // Create FArrayBox for primitive variables
    amrex::Box prims_box(amrex::IntVect(0,0,0), amrex::IntVect(1,1,1));
    amrex::FArrayBox prims_fab(prims_box, ncomp_prims);
    auto prims = prims_fab.array();

    // Initialize primitive variables
    prims(i, j, k, QRHO) = 1.2;      // kg/m^3
    prims(i, j, k, QU) = 100.0;      // m/s
    prims(i, j, k, QV) = 50.0;       // m/s
#if (AMREX_SPACEDIM == 3)
    prims(i, j, k, QW) = 25.0;       // m/s
#endif
    prims(i, j, k, QPRES) = 101325.0; // Pa
    prims(i, j, k, QTEMP) = 300.0;    // K

    // Species mass fractions
    prims(i, j, k, QY1) = 0.767;     // O2 mass fraction
    prims(i, j, k, QY2) = 0.233;     // N2 mass fraction

    // Create dx array
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx;
    dx[0] = 0.01;  // 1 cm
    dx[1] = 0.01;  // 1 cm
#if (AMREX_SPACEDIM == 3)
    dx[2] = 0.01;  // 1 cm
#endif

    // Create thermo object
    pyro::pyro<double> thermo;

    // Call the function
    amrex::Real dt = cns_computedt(i, j, k, prims, dx, thermo);

    // Verify results
    REQUIRE(dt > 0.0);
    REQUIRE(dt < 1.0e10);  // Should be less than the initial large value

    // The time step should be limited by the CFL condition
    // dt = dx / (|u| + c), where c is the speed of sound
    // For a reasonable speed of sound (~340 m/s for air),
    // dt should be approximately dx / (|u| + c)
    double expected_dt = dx[0] / (100.0 + 340.0);  // Rough estimate
    REQUIRE(dt < expected_dt * 2.0);  // Allow some tolerance

    std::cout << "  ✓ cns_computedt test passed" << std::endl;
    std::cout << "    Computed dt: " << dt << std::endl;
    std::cout << "    Expected order of magnitude: " << expected_dt << std::endl;
}

// Round-trip test: cons -> prims -> cons
TEST_CASE_METHOD(AMReXFixture, "roundtrip_conversion", "[thermo]") {
    std::cout << "Testing round-trip conversion (cons -> prims -> cons)..." << std::endl;

    // Create test data
    const int i = 0, j = 0, k = 0;
    const int ncomp_cons = NCONS;
    const int ncomp_prims = NPRIM;

    // Create FArrayBoxes
    amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(1,1,1));
    amrex::FArrayBox cons_original_fab(box, ncomp_cons);
    amrex::FArrayBox cons_converted_fab(box, ncomp_cons);
    amrex::FArrayBox prims_fab(box, ncomp_prims);

    auto cons_original = cons_original_fab.array();
    auto cons_converted = cons_converted_fab.array();
    auto prims = prims_fab.array();

    // Initialize original conservative variables
    cons_original(i, j, k, URHO) = 1.2;
    cons_original(i, j, k, UMX) = 1.2 * 100.0;
    cons_original(i, j, k, UMY) = 1.2 * 50.0;
#if (AMREX_SPACEDIM == 3)
    cons_original(i, j, k, UMZ) = 1.2 * 25.0;
#endif
    cons_original(i, j, k, UEINT) = 1.2 * 200000.0 + 0.5 * 1.2 * (100.0*100.0 + 50.0*50.0 + 25.0*25.0);
    cons_original(i, j, k, UY1) = 1.2 * 0.767;
    cons_original(i, j, k, UY2) = 1.2 * 0.233;

    // Initialize primitive variables with initial guess
    prims(i, j, k, QTEMP) = 300.0;

    // Create thermo object
    pyro::pyro<double> thermo;

    // Convert cons -> prims
    cns_ctoprim(i, j, k, cons_original, prims, thermo);

    // Convert prims -> cons
    cns_primtoc(i, j, k, prims, cons_converted, thermo);

    // Compare original and converted conservative variables
    double tolerance = 1e-6;  // Allow some numerical error for round-trip

    REQUIRE(std::abs(cons_original(i, j, k, URHO) - cons_converted(i, j, k, URHO)) < tolerance);
    REQUIRE(std::abs(cons_original(i, j, k, UMX) - cons_converted(i, j, k, UMX)) < tolerance);
    REQUIRE(std::abs(cons_original(i, j, k, UMY) - cons_converted(i, j, k, UMY)) < tolerance);
#if (AMREX_SPACEDIM == 3)
    REQUIRE(std::abs(cons_original(i, j, k, UMZ) - cons_converted(i, j, k, UMZ)) < tolerance);
#endif
    REQUIRE(std::abs(cons_original(i, j, k, UEINT) - cons_converted(i, j, k, UEINT)) < tolerance);
    REQUIRE(std::abs(cons_original(i, j, k, UY1) - cons_converted(i, j, k, UY1)) < tolerance);
    REQUIRE(std::abs(cons_original(i, j, k, UY2) - cons_converted(i, j, k, UY2)) < tolerance);

    std::cout << "  ✓ Round-trip test passed" << std::endl;
    std::cout << "    Original density: " << cons_original(i, j, k, URHO) << std::endl;
    std::cout << "    Converted density: " << cons_converted(i, j, k, URHO) << std::endl;
    std::cout << "    Difference: " << std::abs(cons_original(i, j, k, URHO) - cons_converted(i, j, k, URHO)) << std::endl;
}

// Test edge cases
TEST_CASE_METHOD(AMReXFixture, "edge_cases", "[thermo]") {
    std::cout << "Testing edge cases..." << std::endl;

    // Test with very small density
    const int i = 0, j = 0, k = 0;
    const int ncomp_cons = NCONS;
    const int ncomp_prims = NPRIM;

    amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(1,1,1));
    amrex::FArrayBox cons_fab(box, ncomp_cons);
    amrex::FArrayBox prims_fab(box, ncomp_prims);

    auto cons = cons_fab.array();
    auto prims = prims_fab.array();

    // Test with very small density
    cons(i, j, k, URHO) = 1e-10;
    cons(i, j, k, UMX) = 1e-10 * 100.0;
    cons(i, j, k, UMY) = 1e-10 * 50.0;
#if (AMREX_SPACEDIM == 3)
    cons(i, j, k, UMZ) = 1e-10 * 25.0;
#endif
    cons(i, j, k, UEINT) = 1e-10 * 200000.0;
    cons(i, j, k, UY1) = 1e-10 * 0.767;
    cons(i, j, k, UY2) = 1e-10 * 0.233;

    prims(i, j, k, QTEMP) = 300.0;

    pyro::pyro<double> thermo;

    // This should handle the small density case gracefully
    cns_ctoprim(i, j, k, cons, prims, thermo);

    // Check that the function doesn't crash and produces reasonable results
    REQUIRE(prims(i, j, k, QRHO) >= 0.0);
    REQUIRE(prims(i, j, k, QPRES) >= 0.0);
    REQUIRE(prims(i, j, k, QTEMP) >= 0.0);

    std::cout << "  ✓ Edge case test passed" << std::endl;
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
