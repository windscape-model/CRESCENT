#include "CRESCENT_NestedSetup.H"

// Primary nested-L1 initialization implementation.

#include <AMReX_Print.H>

#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedSetupDerivedFields.H"
#include "CRESCENT_NestedSetupStorage.H"
#include "CRESCENT_NestedSetupWaveState.H"

namespace crescent::nested::setup {

void initializeNestedL1Level(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    const auto core = level_context.coreAccess();
    const auto nested = level_context.nestedAccess();
    AMREX_ALWAYS_ASSERT(crescent::nested::isActiveNestedL1Level(
        nested.config,
        lev));

    storage::initializeDomainsAndStorage(level_context);
    wave_state::initializeFromCoreParent(level_context);
    derived_fields::initializeFromState(level_context);

    if (core.runtime.verbosity >= 1) {
        amrex::Print() << "=== Nested AMR buffer zone setup complete ===" << std::endl;
    }
}

} // namespace crescent::nested::setup
