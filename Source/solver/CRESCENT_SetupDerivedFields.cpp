#include "CRESCENT_SetupDerivedFields.H"

#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_Physics.H"

namespace crescent::solver::setup::derived_fields {

void initializeCoreDerivedFields(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    AMREX_ALWAYS_ASSERT(lev == 0);
    const auto core = level_context.coreAccess();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            core,
            level_context.levelMesh().domain());
    auto derivative_bundle = core.solution.derivativeBundleAt(lev);
    auto& solution = core.solution;

    crescent::dno::computeSolutionStateDerivativeBundle(
        level_wave_operator,
        core,
        derivative_bundle);

    crescent::physics::computeHorizontalVelocities(
        solution.u[lev],
        solution.v[lev],
        derivative_bundle,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.spectral_workspace.dealiasing_scratch[level_wave_operator.lev]);
}

} // namespace crescent::solver::setup::derived_fields
