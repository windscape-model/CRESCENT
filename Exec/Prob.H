#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#include "ProbBase.H"

class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* problo, const amrex::Real* probhi);

    void init_data(const amrex::Box &bx,
                const amrex::GeometryData &geomdata,
                amrex::Array4<amrex::Real> const &eta,
                amrex::Array4<amrex::Real> const &phi) override;

};

#endif
