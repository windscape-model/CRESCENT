#include "CRESCENT_Init.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_Spectral.H"

#include <AMReX_GpuComplex.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cmath>

using namespace amrex;
using namespace crescent::constants;

namespace {

struct JonswapPeakData
{
    Real alpha;
    Real peak_frequency;
    Real peak_wavelength;
};

JonswapPeakData computeJonswapPeakData(const Real wind_speed_10m, const Real fetch)
{
    JonswapPeakData peak_data;
    peak_data.alpha = 0.076 * std::pow((wind_speed_10m*wind_speed_10m) / fetch / GRAVITY, 0.22);
    peak_data.peak_frequency = 22.0 * std::pow((GRAVITY*GRAVITY / wind_speed_10m / fetch),(1.0/3.0));
    peak_data.peak_wavelength = PI2 * GRAVITY / (peak_data.peak_frequency*peak_data.peak_frequency);
    return peak_data;
}

void initializeJonswapSpectrum(cMultiFab& eta1,
                               cMultiFab& phi1,
                               const int nx,
                               const int ny,
                               const Real dk_x,
                               const Real dk_y,
                               const Real inverse_gravity,
                               const Real alpha,
                               const Real peak_frequency,
                               const Real initial_angle,
                               const Real angle_spread,
                               const int directional_exponent)
{
    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    for (MFIter mfi(eta1, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();
        bx.setRange(0, 0, nx/2 + 1);
        bx.setRange(1, 0, ny/2);

        auto eta1_arr = eta1.array(mfi);
        auto phi1_arr = phi1.array(mfi);

        const Real dk_x_ = dk_x;
        const Real dk_y_ = dk_y;

        ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int z, const RandomEngine& engine) noexcept
        {
            if (i > 0 || j > 0) {
                Real k_x = i * dk_x_;
                Real k_y = j * dk_y_;
                Real k_mag = std::sqrt(k_x*k_x + k_y*k_y);

                Real theta_pos = amrex::Random(engine) * PI2;
                Real theta_neg = amrex::Random(engine) * PI2;

                Real omega = sqrt(k_mag / inverse_gravity);
                Real sigma;
                if (omega < peak_frequency) {
                    sigma = 0.07;
                } else {
                    sigma = 0.09;
                }

                Real peak_enhancement = std::exp(-std::pow((omega - peak_frequency),2)/ 2.0 / std::pow(sigma,2) / std::pow(peak_frequency,2));
                Real jonswap_base = alpha / (inverse_gravity*inverse_gravity) / std::pow(omega,5) * std::exp(-1.25*std::pow((omega/peak_frequency),-4.0)) * std::pow(3.3, peak_enhancement);
                jonswap_base = jonswap_base / (inverse_gravity*inverse_gravity) / 2.0 / std::pow(omega, 3);

                Real directional_factor_pos, directional_factor_neg;

                if (angle_spread < 1e-12) {
                    directional_factor_pos = (j == 0) ? 1.0 : 0.0;
                    directional_factor_neg = (j == 0) ? 1.0 : 0.0;
                } else {
                    Real wave_angle = (j > 0) ? std::atan2(k_y, k_x) : 0.0;
                    if (wave_angle < 0.0) wave_angle += PI2;
                    Real angle_diff = wave_angle - initial_angle;
                    while (angle_diff > PI) angle_diff -= PI2;
                    while (angle_diff < -PI) angle_diff += PI2;

                    Real abs_angle_diff = std::abs(angle_diff);

                    if (abs_angle_diff <= 0.5*angle_spread) {
                        Real directional_factor = 2.0/angle_spread * std::pow(std::cos(PI2/2.0*abs_angle_diff/angle_spread), directional_exponent);
                        directional_factor_pos = directional_factor;
                        directional_factor_neg = directional_factor;
                    } else {
                        directional_factor_pos = 0.0;
                        directional_factor_neg = 0.0;
                    }
                }

                Real spectrum_value_pos = jonswap_base * directional_factor_pos;
                Real spectrum_value_neg = jonswap_base * directional_factor_neg;

                Real wave_amplitude_pos = sqrt(2.0 * spectrum_value_pos * dk_x_ * dk_y_);
                Real wave_amplitude_neg = sqrt(2.0 * spectrum_value_neg * dk_x_ * dk_y_);

                wave_amplitude_pos *= 0.25;
                wave_amplitude_neg *= 0.25;
                if (i == 0 || j == 0) {
                    wave_amplitude_pos *= sqrt(2.0);
                    wave_amplitude_neg *= sqrt(2.0);
                }

                Real phi_scale_pos = wave_amplitude_pos / inverse_gravity / omega;
                Real phi_scale_neg = wave_amplitude_neg / inverse_gravity / omega;

                Real cos_theta_pos = std::cos(theta_pos);
                Real sin_theta_pos = std::sin(theta_pos);
                Real cos_theta_neg = std::cos(theta_neg);
                Real sin_theta_neg = std::sin(theta_neg);

                GpuComplex<Real> eta_pos = GpuComplex<Real>(
                    wave_amplitude_pos * cos_theta_pos,
                    wave_amplitude_pos * sin_theta_pos);
                GpuComplex<Real> phi_pos = GpuComplex<Real>(
                    +phi_scale_pos * sin_theta_pos,
                    -phi_scale_pos * cos_theta_pos);

                GpuComplex<Real> eta_neg = GpuComplex<Real>(
                    wave_amplitude_neg * cos_theta_neg,
                    wave_amplitude_neg * sin_theta_neg);
                GpuComplex<Real> phi_neg = GpuComplex<Real>(
                    +phi_scale_neg * sin_theta_neg,
                    -phi_scale_neg * cos_theta_neg);

                eta1_arr(i, j, z) = eta_pos;
                phi1_arr(i, j, z) = phi_pos;

                if (j > 0) {
                    eta1_arr(i, ny - j, z) = eta_neg;
                    phi1_arr(i, ny - j, z) = phi_neg;
                }
            }
        });
    }
}

} // namespace

namespace crescent::init {

void generateJonswapWaves(MultiFab& eta_saved,
                          MultiFab& phi_saved,
                          std::unique_ptr<FFT::R2C<Real>>& fft,
                          const GeometryData& geom,
                          const BoxArray& cba_root,
                          const DistributionMapping& cdmap_root,
                          const Real dk_x,
                          const Real dk_y,
                          const Real wind_speed_10m,
                          const Real fetch,
                          const Real inverse_gravity,
                          const Real initial_angle,
                          const Real angle_spread,
                          const int directional_exponent,
                          Real& alpha,
                          Real& peak_frequency,
                          Real& peak_wavelength)
{
    const JonswapPeakData peak_data = computeJonswapPeakData(wind_speed_10m, fetch);
    alpha = peak_data.alpha;
    peak_frequency = peak_data.peak_frequency;
    peak_wavelength = peak_data.peak_wavelength;

    auto const &[cba, cdm] = fft->getSpectralDataLayout();
    cMultiFab eta1(cba_root, cdmap_root, 1, 0);
    cMultiFab phi1(cba_root, cdmap_root, 1, 0);

    eta1.setVal(0.0);
    phi1.setVal(0.0);

    cMultiFab eta0_c(cba, cdm, 1, 0);
    cMultiFab phi0_c(cba, cdm, 1, 0);

    const int nx = geom.Domain().length(0);
    const int ny = geom.Domain().length(1);

    initializeJonswapSpectrum(
        eta1, phi1, nx, ny, dk_x, dk_y, inverse_gravity,
        alpha,
        peak_frequency,
        initial_angle,
        angle_spread,
        directional_exponent);

    ParallelDescriptor::Barrier();

    ParallelCopy(eta0_c, eta1, 0, 0, 1);
    ParallelCopy(phi0_c, phi1, 0, 0, 1);

    crescent::spectral::backwardFFT(fft, geom, eta0_c, eta_saved);
    crescent::spectral::backwardFFT(fft, geom, phi0_c, phi_saved);

    // Init-time local scratch; hot-path callers use the persistent
    // SpectralWorkspace dealiasing buffer instead.
    amrex::MultiFab dealiasing_scratch(
        eta_saved.boxArray(),
        eta_saved.DistributionMap(),
        1,
        eta_saved.nGrow());
    crescent::spectral::dealiasing(eta_saved, fft, geom, dealiasing_scratch);
    crescent::spectral::dealiasing(phi_saved, fft, geom, dealiasing_scratch);
}

} // namespace crescent::init
