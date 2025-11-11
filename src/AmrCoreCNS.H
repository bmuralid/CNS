#ifndef AmrCoreCNS_H_
#define AmrCoreCNS_H_

#include <string>
#include <limits>
#include <memory>

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>
#include <AMReX_EBFabFactory.H>
#include <pyro.H>

class AmrCoreCNS
    : public amrex::AmrCore
{
public:

    AmrCoreCNS ();
    virtual ~AmrCoreCNS();

    void InitData ();
    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    void Evolve ();

    static inline pyro::pyro<double> thermo;
    void Cons2Prims(int opt);

    void Prims2Cons();
    void Advance (amrex::Real time, amrex::Real dt);

    void Sensor (amrex::Real time, amrex::Real dt);

    void AdvanceSingleStage (amrex::Real time, amrex::Real dt, int istage);

    amrex::Real ComputeTimeStep (amrex::Real time);

private:

    void ReadParameters();

    void AverageDown(int opt);

    void AverageDownTo (int crse_lev, int opt);


    void FillPatch (int lev, int opt, amrex::Real time, amrex::MultiFab& mf, int  icomp, int ncomp);

    void FillCoarsePatch (int lev, int opt, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    void GetData (int lev, int opt,
                  amrex::Real time,
                  amrex::Vector<amrex::MultiFab*>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    // void ComputeTimeStep ();

    void TimeStep (amrex::Real time);


    // private data members
    amrex::Vector<int> istep;
    amrex::Vector<int> nsubsteps;

    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    amrex::Vector<amrex::MultiFab> qcons_new;
    amrex::Vector<amrex::MultiFab> qcons_tmp;
    amrex::Vector<amrex::MultiFab> qcons_old;
    amrex::Vector<amrex::MultiFab> qprims;
    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> dq;
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory> > ebfactory;

    amrex::Vector<amrex::BCRec> bcs;
    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg;
    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // advective cfl number - dt = cfl*dx/umax
    amrex::Real cfl = 0.7;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 1;

    // do we subcycle in time?
    int do_subcycle = 1;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;
};

#endif
