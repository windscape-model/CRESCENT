#pragma once

// Primary solver declaration header.

#include <memory>

#include <AMReX_AmrCore.H>
#include "CRESCENT_ExternalExchange.H"

namespace crescent {
struct SolverState;
}

namespace crescent::shared {
struct MeshHierarchyView;
class SolverContext;
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
    [[nodiscard]] crescent::SolverState& solverState() noexcept;
    [[nodiscard]] const crescent::SolverState& solverState() const noexcept;
    [[nodiscard]] crescent::shared::MeshHierarchyView meshHierarchy() const;
    [[nodiscard]] crescent::shared::SolverContext solverContext();

    std::unique_ptr<crescent::SolverState> solver_state;
};
