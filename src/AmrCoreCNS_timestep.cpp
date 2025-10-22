
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Interp_C.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <cmath>

#include <AmrCoreCNS.H>
#include <Kernels.H>

using namespace amrex;

amrex::Real
AmrCoreCNS::ComputeTimeStep (Real time)
{
    pyro::pyro<double> const lpyro  = thermo;

    Real dt_min = std::numeric_limits<Real>::max();

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (int lev = 0; lev <= finest_level; lev++) {
        MultiFab& mfprims =  qprims[lev];

        const auto dx = Geom(lev).CellSizeArray();
        for (MFIter mfi(mfprims, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real> pfab = mfprims[mfi].array();

            const Box& bx = mfi.tilebox();
            reduce_op.eval(bx, reduce_data, [=]
            AMREX_GPU_DEVICE (int i, int j, int k)  -> ReduceTuple
            {
                Real dt_loc = cns_computedt(i, j, k, pfab, dx, lpyro);
                return {dt_loc};
            });

        }
    }

    ReduceTuple hv = reduce_data.value();
    dt_min = amrex::min(dt_min, amrex::get<0>(hv));

    dt_min *= cfl;
    ParallelDescriptor::ReduceRealMin(dt_min);
    return dt_min;
}

