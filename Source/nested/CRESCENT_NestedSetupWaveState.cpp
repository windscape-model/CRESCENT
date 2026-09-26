#include "CRESCENT_NestedSetupWaveState.H"

#include <AMReX_Print.H>

#include "CRESCENT_Interpolation.H"

namespace crescent::nested::setup::wave_state {

void initializeFromCoreParent(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    AMREX_ALWAYS_ASSERT(lev > 0);
    const auto core = level_context.coreAccess();
    const auto nested = level_context.nestedAccess();
    const int verbosity = level_context.runtime().verbosity;
    if (verbosity >= 1) {
        amrex::Print() << "\n=== Initializing L1 from parent L0 ===" << std::endl;
    }

    const auto core_mesh = level_context.coreLevelMesh();
    const auto level_mesh = level_context.levelMesh();
    auto parent_state = core.solution.freeSurfaceStateAt(0);
    auto child_state = core.solution.freeSurfaceStateAt(lev);
    auto& coarse_interp_patch =
        nested.runtime.post_rk4_state_scratch[lev].coarseInterpPatch();

    parent_state.eta.FillBoundary(core_mesh.geom.periodicity());
    parent_state.phi.FillBoundary(core_mesh.geom.periodicity());

    crescent::interpolation::interpolateParentToChild(
        crescent::shared::asConstFreeSurfaceState(parent_state),
        child_state,
        core_mesh.geom,
        level_mesh.geom,
        level_context.refineRatio(),
        coarse_interp_patch,
        nested.config.state_interpolation_order);

    core.solution.saveSurfaceState(lev);
}

} // namespace crescent::nested::setup::wave_state
