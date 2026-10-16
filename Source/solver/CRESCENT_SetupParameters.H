#pragma once

namespace crescent {
struct SolverState;
}

namespace crescent::solver::setup {

// Reads and validates every input of the "crescent" ParmParse table (and
// max_step, start_time and stop_time) into the solver state.
void readParameters(crescent::SolverState& state, int max_level);

} // namespace crescent::solver::setup
