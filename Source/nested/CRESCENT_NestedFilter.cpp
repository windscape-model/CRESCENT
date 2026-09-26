#include "CRESCENT_NestedFilter.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_SpectralFilterExecution.H"
#include "CRESCENT_SpectralGridMetrics.H"
#include "CRESCENT_SpectralFilterKernels.H"

// Active nested-L1 filter implementation.

namespace crescent::nested::filter {

namespace FFT = amrex::FFT;

using amrex::GeometryData;
using amrex::GpuComplex;
using amrex::MultiFab;
using amrex::Real;

void applyHouLiFilter_L1(
    MultiFab& eta,
    MultiFab& phi,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geomdata,
    Real dk_x,
    Real dk_y,
    Real k_peak,
    Real beta1,
    Real beta2)
{
    // Get grid dimensions
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(eta, geomdata);
    const int nx = spectral_grid.nx;
    const int ny = spectral_grid.ny;

    // CRITICAL FIX: Compute k_nyquist from dk_x (wavenumber spacing), not from geomdata!
    // For L1, geomdata gives L0 domain size, which would give wrong dx and wrong k_nyquist.
    // k_nyquist = (nx/2) * dk_x, since the highest mode is nx/2 with spacing dk_x.
    const Real k_nyquist = spectral_grid.k_nyquist_x;
    const Real k_dealiasing = (2.0 / 3.0) * k_nyquist;
    const Real normalization_factor = Real(nx * ny);

    const crescent::spectral::filter_kernels::CombinedSpectralFilterParameters
        filter_params{
            ny,
            dk_x,
            dk_y,
            normalization_factor,
            k_nyquist,
            k_dealiasing,
            k_peak,
            beta1,
            beta2,
            crescent::constants::HOULIGHER_FILTER_ALPHA_AGGRESSIVE,
            static_cast<int>(crescent::constants::HOULIGHER_FILTER_POWER),
            true,
            true
        };

    const auto spectral_filter =
        [=] AMREX_GPU_DEVICE(int i, int j, int k, GpuComplex<Real>& fft_c) noexcept
        {
            crescent::spectral::filter_kernels::
                applyCombinedHouLiDealiasLowPassZeroMean(
                    i,
                    j,
                    fft_c,
                    filter_params);
        };

    crescent::spectral::applyInPlaceSpectralFilter(eta, fft, spectral_filter);
    crescent::spectral::applyInPlaceSpectralFilter(phi, fft, spectral_filter);
}

} // namespace crescent::nested::filter
