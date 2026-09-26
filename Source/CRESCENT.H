#pragma once

// Primary solver declaration header.

#include <memory>

#include <AMReX_AmrCore.H>
#include "CRESCENT_ExternalExchange.H"

namespace crescent::ownership {
struct SolverOwnership;
}

class CRESCENT : public amrex::AmrCore,
               public crescent::external_exchange::SurfaceExchangeProvider
{
public:
    CRESCENT();
    ~CRESCENT() override;

    CRESCENT(CRESCENT&&) noexcept = delete;
    CRESCENT& operator=(CRESCENT&& other) noexcept = delete;
    CRESCENT(const CRESCENT& other) = delete;
    CRESCENT& operator=(const CRESCENT& other) = delete;

    void Init();
    void Run();

    [[nodiscard]] int finestSurfaceExchangeLevel() const override;
    [[nodiscard]] crescent::external_exchange::FreeSurfaceExchangeView surfaceExchangeView(int lev) override;
    [[nodiscard]] crescent::external_exchange::ConstFreeSurfaceExchangeView surfaceExchangeView(int lev) const override;

    // ------------------------------------
    // Members to override from AmrCore:
    // ------------------------------------
    // Tag cells for refinement
    void ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time, int ngrow) override;

    void MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray &ba,
                                const amrex::DistributionMapping &dm) override;
    void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray &ba,
                     const amrex::DistributionMapping &dm) override;
    void ClearLevel(int lev) override;
    void MakeNewLevelFromScratch(int lev, amrex::Real time, const amrex::BoxArray &ba,
                                 const amrex::DistributionMapping &dm) override;
    // ------------------------------------

private:
    [[nodiscard]] crescent::ownership::SolverOwnership& solverOwnership() noexcept;
    [[nodiscard]] const crescent::ownership::SolverOwnership& solverOwnership() const noexcept;

    std::unique_ptr<crescent::ownership::SolverOwnership> solver_ownership;
};
