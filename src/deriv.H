#ifndef THERMO_H
#define THERMO_H

#include <AMReX_Array4.H>
#include <pyro.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_deri (int i, int j, int k,
        amrex::Array4<amrex::Real> const& cons,
        amrex::Array4<amrex::Real> const& dcons,
        pyro::pyro<double> const& thermo)
{


    dcons(i, j, k, 0)  = cons(i+1, j, k, 0) - cons(i-1, j, k, 0) ;




}
