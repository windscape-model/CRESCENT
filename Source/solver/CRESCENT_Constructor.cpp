#include "CRESCENT.H"
#include "CRESCENT_SetupLifecycle.H"
#include "CRESCENT_SetupParameters.H"
#include "CRESCENT_SolverOwnership.H"

// Primary implementation file.

CRESCENT::CRESCENT()
    : solver_ownership(std::make_unique<crescent::ownership::SolverOwnership>())
{
    const auto construction_context =
        crescent::solver::setup::makeConstructionContext(
            max_level,
            finest_level,
            solverOwnership().coreAccess(),
            solverOwnership().nestedAccess(),
            solverOwnership().amrAccess(),
            geom[0]);
    crescent::solver::setup::initializeLevelTimeState(
        construction_context);
    crescent::solver::setup::readParameters(construction_context);
    crescent::solver::setup::initializeConfiguredConstructionState(
        construction_context);
}

crescent::ownership::SolverOwnership& CRESCENT::solverOwnership() noexcept
{
    return *solver_ownership;
}

const crescent::ownership::SolverOwnership& CRESCENT::solverOwnership() const noexcept
{
    return *solver_ownership;
}
