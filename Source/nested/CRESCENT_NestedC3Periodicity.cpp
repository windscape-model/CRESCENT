#include "CRESCENT_NestedC3Periodicity.H"
#include "CRESCENT_NestedC3LinearSolvers.H"

// Primary nested C3-repair implementation.
#include "CRESCENT_NestedDomains.H"
#include <AMReX_Gpu.H>
#include <cmath>

namespace crescent::nested::c3 {

using namespace amrex;

// ============================================================================
// C∞ smooth step function (from user's plan)
// ============================================================================

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real smooth_phi(Real s) {
    if (s <= 0.0) return 0.0;
    return std::exp(-1.0 / s);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real smooth_step(Real s) {
    // S(s) in [0,1], C∞, flat at 0 and 1
    Real a = smooth_phi(s);
    Real b = smooth_phi(1.0 - s);
    return (a + b > 0.0) ? a / (a + b) : 0.0;
}

// ============================================================================
// Factorial helper
// ============================================================================

inline Real factorial(int n) {
    Real f = 1.0;
    for (int i = 2; i <= n; ++i) f *= i;
    return f;
}

// Legacy polynomial basis kept for applyC3SeamRepair.
inline Real p_r(Real xi, int r) {
    Real one_minus_xi = 1.0 - xi;
    Real v = 1.0;
    for (int m = 0; m < r + 1; ++m) v *= one_minus_xi;
    return v / factorial(r + 1);
}

// ============================================================================
// 4th-order one-sided derivative stencils
// ============================================================================

void compute_left_derivatives_1d(const Real* f, Real dx, Real* deriv) {
    // Value
    deriv[0] = f[0];

    // 1st derivative
    deriv[1] = (-25.0*f[0] + 48.0*f[1] - 36.0*f[2] + 16.0*f[3] - 3.0*f[4]) / (12.0*dx);

    // 2nd derivative
    deriv[2] = (35.0*f[0] - 104.0*f[1] + 114.0*f[2] - 56.0*f[3] + 11.0*f[4]) / (12.0*dx*dx);

    // 3rd derivative
    deriv[3] = (-5.0*f[0] + 18.0*f[1] - 24.0*f[2] + 14.0*f[3] - 3.0*f[4]) / (2.0*dx*dx*dx);
}

void compute_right_derivatives_1d(const Real* f, int n, Real dx, Real* deriv) {
    // Mirror indexing
    Real g[5];
    for (int m = 0; m < 5; ++m) {
        g[m] = f[n-1-m];
    }

    // Value
    deriv[0] = g[0];

    // 1st derivative (sign flip)
    deriv[1] = (25.0*g[0] - 48.0*g[1] + 36.0*g[2] - 16.0*g[3] + 3.0*g[4]) / (12.0*dx);

    // 2nd derivative
    deriv[2] = (35.0*g[0] - 104.0*g[1] + 114.0*g[2] - 56.0*g[3] + 11.0*g[4]) / (12.0*dx*dx);

    // 3rd derivative (sign flip)
    deriv[3] = (5.0*g[0] - 18.0*g[1] + 24.0*g[2] - 14.0*g[3] + 3.0*g[4]) / (2.0*dx*dx*dx);
}

// ============================================================================
// Shared linear-system and modal solver kernels live in
// CRESCENT_NestedC3LinearSolvers.cpp.

// ============================================================================
// Main C³-repair function with ρ-based correction
// ============================================================================

void applyC3SeamRepair(MultiFab& field,
                       const Geometry& geom,
                       const BoxArray& ba,
                       int buffer_width)
{
    const int k = 8;  // Number of basis functions (overdetermined)

    // Get GLOBAL L1 box dimensions from BoxArray minimal box
    // For multi-box layouts, we need the encompassing box
    const Box global_box = ba.minimalBox();
    int nx = global_box.length(0);
    int ny = global_box.length(1);

    // Cell size
    Real dx = geom.CellSize(0);
    Real dy = geom.CellSize(1);

    // L1 box physical size
    Real Lx = nx * dx;
    Real Ly = ny * dy;

    // Safety checks
    AMREX_ALWAYS_ASSERT(nx >= 5 && ny >= 5);

    int bx = buffer_width;
    int by = buffer_width;
    // Minimum buffer: 5 points for stencil + a few for smooth transition
    AMREX_ALWAYS_ASSERT(bx >= 8 && by >= 8);

    // ========================================================================
    // STEP 1: Build ρ_x and ρ_y masks (0 in ROI, ramp to 1 at seams)
    // ========================================================================

    Vector<Real> rho_x(nx, 0.0);
    for (int i = 0; i < nx; ++i) {
        if (i < bx) {
            // Left buffer: s=1 at seam (i=0), s=0 at ROI boundary (i=bx)
            Real s = Real(bx - i) / Real(bx);
            rho_x[i] = smooth_step(s);
        } else if (i > nx - 1 - bx) {
            // Right buffer
            int i0 = nx - bx - 1;  // Inner interface
            Real s = Real(i - i0) / Real(bx);
            rho_x[i] = smooth_step(s);
        } else {
            rho_x[i] = 0.0;  // ROI: exactly 0
        }
    }

    Vector<Real> rho_y(ny, 0.0);
    for (int j = 0; j < ny; ++j) {
        if (j < by) {
            Real s = Real(by - j) / Real(by);
            rho_y[j] = smooth_step(s);
        } else if (j > ny - 1 - by) {
            int j0 = ny - by - 1;
            Real s = Real(j - j0) / Real(by);
            rho_y[j] = smooth_step(s);
        } else {
            rho_y[j] = 0.0;
        }
    }

    // ========================================================================
    // STEP 2: Build basis Q_x[r][i] = rho_x[i] * p_r(xi_i)
    // ========================================================================

    std::vector<std::vector<Real>> Qx(k, std::vector<Real>(nx, 0.0));
    for (int r = 0; r < k; ++r) {
        for (int i = 0; i < nx; ++i) {
            Real xi = (Real(i) + 0.5) / Real(nx);  // Normalized coordinate
            Qx[r][i] = rho_x[i] * p_r(xi, r);
        }
    }

    std::vector<std::vector<Real>> Qy(k, std::vector<Real>(ny, 0.0));
    for (int r = 0; r < k; ++r) {
        for (int j = 0; j < ny; ++j) {
            Real xi = (Real(j) + 0.5) / Real(ny);
            Qy[r][j] = rho_y[j] * p_r(xi, r);
        }
    }

    // ========================================================================
    // STEP 3: Build jump matrices B_x and B_y
    // ========================================================================

    std::vector<std::vector<Real>> Bx(4, std::vector<Real>(k, 0.0));
    for (int r = 0; r < k; ++r) {
        Real DL[4], DR[4];
        compute_left_derivatives_1d(Qx[r].data(), dx, DL);
        compute_right_derivatives_1d(Qx[r].data(), nx, dx, DR);
        for (int n = 0; n < 4; ++n) {
            Bx[n][r] = DR[n] - DL[n];
        }
    }

    std::vector<std::vector<Real>> By(4, std::vector<Real>(k, 0.0));
    for (int r = 0; r < k; ++r) {
        Real DL[4], DR[4];
        compute_left_derivatives_1d(Qy[r].data(), dy, DL);
        compute_right_derivatives_1d(Qy[r].data(), ny, dy, DR);
        for (int n = 0; n < 4; ++n) {
            By[n][r] = DR[n] - DL[n];
        }
    }

    // ========================================================================
    // STEP 4: Gather multi-box data into a single-box scratch MultiFab
    // ========================================================================
    BoxArray ba_global(global_box);
    DistributionMapping dm_global(ba_global);
    MultiFab field_global(ba_global, dm_global, 1, 0);

    // Gather: Copy from multi-box field to single-box field_global
    field_global.ParallelCopy(field, 0, 0, 1);

    // Now process the single-box field_global using the original logic
    for (MFIter mfi(field_global); mfi.isValid(); ++mfi) {
        const Box& bx_box = mfi.validbox();
        auto const& field_arr = field_global.array(mfi);

        // Verify this is the global box
        AMREX_ALWAYS_ASSERT(bx_box == global_box);

        // Extract 2D array
        std::vector<Real> f_2d(nx * ny);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = j * nx + i;
                f_2d[idx] = field_arr(i + bx_box.smallEnd(0),
                                     j + bx_box.smallEnd(1), 0);
            }
        }

        // X-PASS
        for (int j = 0; j < ny; ++j) {
            // Extract line (NO corner masking - ρ already in basis)
            Vector<Real> u_line(nx);
            for (int i = 0; i < nx; ++i) {
                u_line[i] = f_2d[j * nx + i];
            }

            // Compute seam jump
            Real DL[4], DR[4], Delta[4];
            compute_left_derivatives_1d(u_line.data(), dx, DL);
            compute_right_derivatives_1d(u_line.data(), nx, dx, DR);
            for (int n = 0; n < 4; ++n) {
                Delta[n] = DR[n] - DL[n];
            }

            // Solve least-squares
            Real alpha[8];
            bool success = solve_least_squares(Bx, Delta, k, alpha);
            if (!success) continue;

            // Apply correction (ρ_x already embedded in Q_x basis)
            for (int i = 0; i < nx; ++i) {
                Real delta = 0.0;
                for (int r = 0; r < k; ++r) {
                    delta += alpha[r] * Qx[r][i];
                }
                u_line[i] += delta;
            }

            // Write back
            for (int i = 0; i < nx; ++i) {
                f_2d[j * nx + i] = u_line[i];
            }
        }

        // Y-PASS
        for (int i = 0; i < nx; ++i) {
            // Extract line (NO corner masking - ρ already in basis)
            Vector<Real> u_line(ny);
            for (int j = 0; j < ny; ++j) {
                u_line[j] = f_2d[j * nx + i];
            }

            // Compute seam jump
            Real DL[4], DR[4], Delta[4];
            compute_left_derivatives_1d(u_line.data(), dy, DL);
            compute_right_derivatives_1d(u_line.data(), ny, dy, DR);
            for (int n = 0; n < 4; ++n) {
                Delta[n] = DR[n] - DL[n];
            }

            // Solve least-squares with Tikhonov regularization
            // When Delta=0, the regularization ensures alpha=0 (no correction)
            Real alpha[8];
            bool success = solve_least_squares(By, Delta, k, alpha);
            if (!success) continue;

            // Apply correction (ρ_y already embedded in Q_y basis)
            for (int j = 0; j < ny; ++j) {
                Real delta = 0.0;
                for (int r = 0; r < k; ++r) {
                    delta += alpha[r] * Qy[r][j];
                }
                u_line[j] += delta;
            }

            // Write back
            for (int j = 0; j < ny; ++j) {
                f_2d[j * nx + i] = u_line[j];
            }
        }

        // Additional sweeps
        // X-PASS (2nd)
        for (int j = 0; j < ny; ++j) {
            Vector<Real> u_line(nx);
            for (int i = 0; i < nx; ++i) u_line[i] = f_2d[j * nx + i];

            Real DL[4], DR[4], Delta[4];
            compute_left_derivatives_1d(u_line.data(), dx, DL);
            compute_right_derivatives_1d(u_line.data(), nx, dx, DR);
            for (int n = 0; n < 4; ++n) Delta[n] = DR[n] - DL[n];

            Real alpha[8];
            if (!solve_least_squares(Bx, Delta, k, alpha)) continue;

            for (int i = 0; i < nx; ++i) {
                Real delta = 0.0;
                for (int r = 0; r < k; ++r) delta += alpha[r] * Qx[r][i];
                u_line[i] += delta;
            }

            for (int i = 0; i < nx; ++i) f_2d[j * nx + i] = u_line[i];
        }

        // Y-PASS (2nd)
        for (int i = 0; i < nx; ++i) {
            Vector<Real> u_line(ny);
            for (int j = 0; j < ny; ++j) u_line[j] = f_2d[j * nx + i];

            Real DL[4], DR[4], Delta[4];
            compute_left_derivatives_1d(u_line.data(), dy, DL);
            compute_right_derivatives_1d(u_line.data(), ny, dy, DR);
            for (int n = 0; n < 4; ++n) Delta[n] = DR[n] - DL[n];

            // Tikhonov regularization ensures alpha=0 when Delta=0
            Real alpha[8];
            if (!solve_least_squares(By, Delta, k, alpha)) continue;

            for (int j = 0; j < ny; ++j) {
                Real delta = 0.0;
                for (int r = 0; r < k; ++r) delta += alpha[r] * Qy[r][j];
                u_line[j] += delta;
            }

            for (int j = 0; j < ny; ++j) f_2d[j * nx + i] = u_line[j];
        }

        // Write back to single-box field_global
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                field_arr(i + bx_box.smallEnd(0),
                         j + bx_box.smallEnd(1), 0) = f_2d[j * nx + i];
            }
        }

    }

    // ========================================================================
    // Scatter - Copy from single-box field_global back to multi-box field
    // ========================================================================
    field.ParallelCopy(field_global, 0, 0, 1);
}


// ============================================================================
void applyC3SeamRepair(MultiFab& field,
                       const Geometry& geom,
                       const BoxArray& ba,
                       const crescent::nested::domains::DomainStructure& domains,
                       int buffer_width)
{
    // The domain structure is accepted for interface parity; the repair uses
    // the provided buffer_width directly.
    (void)domains;
    applyC3SeamRepair(field, geom, ba, buffer_width);
}

} // namespace crescent::nested::c3
