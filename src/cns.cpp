#include <cns.H>

using namespace amrex;

void initializeEB (const Geometry&, int, int);

int cns(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        // timer for profiling
        BL_PROFILE("main()");

        // wallclock time
        const auto strt_total = amrex::second();

        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures
        AmrCoreCNS amr_core_adv;

        // initialize AMR data
        auto max_level = amr_core_adv.maxLevel();
        initializeEB(amr_core_adv.Geom(0), 0 , 0);
        amr_core_adv.InitData();

        // advance solution to final time
        /* amr_core_adv.Evolve(); */

        // wallclock time
        auto end_total = amrex::second() - strt_total;

        if (1) {
        /* if (amr_core_adv.Verbose()) { */
            // print wallclock time
            ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
            amrex::Print() << "\nTotal Time: " << end_total << '\n';
        }
    }

    amrex::Finalize();
    return 0;
}
