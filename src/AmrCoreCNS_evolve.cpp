
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

void
AmrCoreCNS:: Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    /* Actual time step loop */
    for(int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        Real dt = ComputeTimeStep(cur_time);
        Advance(cur_time, dt);
        amrex::Print() << "After step " << step+1 << " at time " << cur_time+dt
                           << " dt: " << dt << "\n";

        cur_time += dt;

        for (int lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

        if (max_level > 0 && regrid_int > 0)
        {
            if (istep[0] % regrid_int == 0)
            {
                regrid(0, cur_time);
            }
        }

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
        }
        if (chk_int > 0 && (step+1)% chk_int == 0) {
            WriteCheckpointFile();
        }
        istep[0] += 1;
    }
    // Final plotfile
    WritePlotFile();
}

