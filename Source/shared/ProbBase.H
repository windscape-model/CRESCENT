#ifndef PROBBASE_H_
#define PROBBASE_H_

#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

class ProblemBase
{
public:
    virtual ~ProblemBase() = default;

    virtual void init_data(const amrex::Box &bx,
                           const amrex::GeometryData &geomdata,
                           amrex::Array4<amrex::Real> const &eta,
                           amrex::Array4<amrex::Real> const &phi) {

                           }
};

extern std::unique_ptr<ProblemBase> amrex_probinit(const amrex_real* problo,
                                                   const amrex_real* probhi) AMREX_ATTRIBUTE_WEAK;

#endif
