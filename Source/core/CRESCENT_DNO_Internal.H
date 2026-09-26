#pragma once

#include "CRESCENT_Constants.H"
#include "CRESCENT_DNOWorkspaceViews.H"
#include "CRESCENT_Spectral.H"
#include "CRESCENT_SpectralFilterKernels.H"
#include <AMReX_FFT.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_MultiFab.H>
#include <cmath>
#include <memory>

namespace crescent::dno::detail {

using amrex::Box;
using amrex::cMultiFab;
using amrex::Copy;
using amrex::GeometryData;
using amrex::GpuComplex;
using amrex::make_alias;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Real;
using amrex::TilingIfNotGPU;

namespace FFT = amrex::FFT;

struct PhiExpansionRecurrenceContext
{
    const MultiFab& surface_phi;
    const MultiFab& wavenumber_pow;
    PhiExpansionWorkspaceView workspace;
    int perturbation_order;
    std::unique_ptr<FFT::R2C<Real>>& fft;
    GeometryData geom;
};

inline PhiExpansionRecurrenceContext makePhiExpansionRecurrenceContext(
    const MultiFab& surface_phi,
    const MultiFab& wavenumber_pow,
    PhiExpansionWorkspaceView workspace,
    int perturbation_order,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom)
{
    return {
        surface_phi,
        wavenumber_pow,
        workspace,
        perturbation_order,
        fft,
        geom
    };
}

struct VerticalVelocityAccumulationContext
{
    MultiFab& phi_z;
    VerticalVelocityWorkspaceView workspace;
    int perturbation_order;
    std::unique_ptr<FFT::R2C<Real>>& fft;
    Real dk_x;
    Real dk_y;
};

inline VerticalVelocityAccumulationContext makeVerticalVelocityAccumulationContext(
    MultiFab& phi_z,
    VerticalVelocityWorkspaceView workspace,
    int perturbation_order,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    Real dk_x,
    Real dk_y)
{
    return {
        phi_z,
        workspace,
        perturbation_order,
        fft,
        dk_x,
        dk_y
    };
}

inline void buildCumulativePhiExpansionCoefficients(
    MultiFab& pk,
    const MultiFab& phi_expansion,
    int perturbation_order)
{
    for (int n = 0; n < perturbation_order - 1; ++n) {
        for (MFIter mfi(pk, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto pk_arr = pk.array(mfi);
            auto phi_arr = phi_expansion.const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                pk_arr(i,j,k,n+1) = pk_arr(i,j,k,n) + phi_arr(i,j,k,n+1);
            });
        }
    }
}

inline void initializeCumulativePhiExpansionState(
    MultiFab& pk,
    const MultiFab& phi_expansion,
    int perturbation_order)
{
    Copy(pk, phi_expansion, 0, 0, 1, 0);

    if (perturbation_order > 1) {
        buildCumulativePhiExpansionCoefficients(pk, phi_expansion, perturbation_order);
    }
}

inline void initializePhiExpansionRecurrence(
    MultiFab& phi_expansion,
    MultiFab& expansion_deriv,
    cMultiFab& deriv_spectral,
    const MultiFab& surface_phi)
{
    expansion_deriv.setVal(0.0);
    deriv_spectral.setVal(0.0);
    Copy(phi_expansion, surface_phi, 0, 0, 1, 0);
}

inline void forwardPhiExpansionComponent(
    cMultiFab& spectral_phi_expansion,
    MultiFab& phi_expansion,
    int component,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom)
{
    MultiFab physical_component(phi_expansion, make_alias, component, 1);
    cMultiFab spectral_component(spectral_phi_expansion, make_alias, component, 1);
    crescent::spectral::forwardFFT(
        fft, geom, physical_component, spectral_component);
}

inline void applyTwoThirdsSpectralCutoff(cMultiFab& spectral_data, int nx, int ny)
{
    int cutoff_x = static_cast<int>(nx * crescent::constants::ONE_THIRD_CUTOFF);
    int cutoff_y = static_cast<int>(ny * crescent::constants::ONE_THIRD_CUTOFF);

    for (MFIter mfi(spectral_data, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.fabbox();
        auto spectral_arr = spectral_data.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const int kx =
                crescent::spectral::filter_kernels::spectralKxIndex(i, nx);
            const int ky =
                crescent::spectral::filter_kernels::spectralKyIndex(j, ny);

            if (crescent::spectral::filter_kernels::spectralModeExceedsCutoff(
                    kx, ky, cutoff_x, cutoff_y)) {
                spectral_arr(i,j,k) = GpuComplex<Real>(0.0, 0.0);
            }
        });
    }
}

inline void finalizePhiExpansionWithTwoThirdsCutoff(
    MultiFab& phi_expansion,
    cMultiFab& spectral_phi_expansion,
    int perturbation_order,
    std::unique_ptr<FFT::R2C<Real>>& fft)
{
    const Box& mf_box = phi_expansion.boxArray().minimalBox();
    const int nx = mf_box.length(0);
    const int ny = mf_box.length(1);

    for (int kw = 0; kw < perturbation_order; ++kw)
    {
        cMultiFab spectral_component(spectral_phi_expansion, make_alias, kw, 1);
        applyTwoThirdsSpectralCutoff(spectral_component, nx, ny);

        MultiFab physical_component(phi_expansion, make_alias, kw, 1);
        fft->backward(spectral_component, physical_component);
    }
}

inline void applyWavenumberPowerInPlace(
    const MultiFab& input,
    MultiFab& output,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    Real dk_x,
    Real dk_y,
    int wvn_power)
{
    const Box& mf_box = input.boxArray().minimalBox();
    const int nx = mf_box.length(0);
    const int ny = mf_box.length(1);
    const Real normalization = static_cast<Real>(nx) * static_cast<Real>(ny);
    const Real dk_x_local = dk_x;
    const Real dk_y_local = dk_y;

    fft->forwardThenBackward(input, output, [=] AMREX_GPU_DEVICE (int i, int j, int k, GpuComplex<Real> &fft_x) noexcept
    {
        const Real k_mag =
            crescent::spectral::filter_kernels::spectralMagnitude(
                i, j, ny, dk_x_local, dk_y_local);
        const Real wvn_val = std::pow(k_mag, wvn_power);

        fft_x = fft_x * wvn_val / normalization;
    });
}

inline void computeVerticalVelocityTermFromPkComponent(
    MultiFab& pk,
    int pk_component,
    MultiFab& term_out,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    Real dk_x,
    Real dk_y,
    int wvn_power)
{
    MultiFab pk_component_mf(pk, make_alias, pk_component, 1);
    applyWavenumberPowerInPlace(
        pk_component_mf, term_out, fft, dk_x, dk_y, wvn_power);
}

inline void multiplySpectralExpansionTerm(
    cMultiFab& deriv_spectral,
    cMultiFab& phi_expansion_spectrum,
    const MultiFab& wavenumber_pow,
    int eta_order,
    int phi_component)
{
    for (MFIter mfi(phi_expansion_spectrum, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.fabbox();

        auto phi_spectrum_arr = phi_expansion_spectrum.array(mfi);
        auto wavenumber_arr = wavenumber_pow.const_array(mfi);
        auto dr_c_arr = deriv_spectral.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dr_c_arr(i,j,k) =
                phi_spectrum_arr(i,j,k,phi_component) *
                wavenumber_arr(i,j,k,eta_order);
        });
    }
}

inline void subtractWeightedPhysicalContribution(
    MultiFab& phi_expansion,
    const MultiFab& expansion_deriv,
    const MultiFab& eta_pow,
    int kw,
    int k1)
{
    for (MFIter mfi(expansion_deriv, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto phi_arr = phi_expansion.array(mfi);
        auto dr_arr = expansion_deriv.const_array(mfi);
        auto eta_pow_arr = eta_pow.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi_arr(i,j,k,kw) = phi_arr(i,j,k,kw) - eta_pow_arr(i,j,k,k1) * dr_arr(i,j,k);
        });
    }
}

inline void advancePhiExpansionOrder(
    PhiExpansionRecurrenceContext& context,
    int expansion_order)
{
    auto& workspace = context.workspace;
    auto& phi_expansion = workspace.phi_expansion;

    phi_expansion.setVal(0.0, expansion_order, 1, 0);

    for (int eta_order = 0; eta_order < expansion_order; ++eta_order)
    {
        const int phi_component = expansion_order - eta_order - 1;

        multiplySpectralExpansionTerm(
            workspace.deriv_spectral,
            workspace.phi_expansion_spectrum,
            context.wavenumber_pow,
            eta_order,
            phi_component);

        crescent::spectral::backwardFFT(
            context.fft,
            context.geom,
            workspace.deriv_spectral,
            workspace.expansion_deriv);

        subtractWeightedPhysicalContribution(
            phi_expansion,
            workspace.expansion_deriv,
            workspace.eta_pow,
            expansion_order,
            eta_order);
    }

    forwardPhiExpansionComponent(
        workspace.phi_expansion_spectrum,
        phi_expansion,
        expansion_order,
        context.fft,
        context.geom);
}

inline void executePhiExpansionRecurrence(
    PhiExpansionRecurrenceContext& context)
{
    auto& workspace = context.workspace;

    initializePhiExpansionRecurrence(
        workspace.phi_expansion,
        workspace.expansion_deriv,
        workspace.deriv_spectral,
        context.surface_phi);
    forwardPhiExpansionComponent(
        workspace.phi_expansion_spectrum,
        workspace.phi_expansion,
        0,
        context.fft,
        context.geom);

    if (context.perturbation_order > 1) {
        for (int expansion_order = 1;
             expansion_order < context.perturbation_order;
             ++expansion_order)
        {
            advancePhiExpansionOrder(
                context,
                expansion_order);
        }
    }

    // OPTIMIZED: Apply 2/3 rule filter in spectral space, then IFFT once
    // This replaces: IFFT + dealiasing(FFT + filter + IFFT) with just: filter + IFFT
    // Reduces 9 FFTs (3 IFFT + 6 from dealiasing) to 3 IFFTs
    finalizePhiExpansionWithTwoThirdsCutoff(
        workspace.phi_expansion,
        workspace.phi_expansion_spectrum,
        context.perturbation_order,
        context.fft);
}

inline void executePhiExpansionRecurrence(
    const MultiFab& surface_phi,
    const MultiFab& wavenumber_pow,
    PhiExpansionWorkspaceView workspace,
    int perturbation_order,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom)
{
    auto context = makePhiExpansionRecurrenceContext(
        surface_phi,
        wavenumber_pow,
        workspace,
        perturbation_order,
        fft,
        geom);
    executePhiExpansionRecurrence(context);
}

inline void accumulateWeightedIntoField(
    MultiFab& output,
    const MultiFab& weights,
    const MultiFab& input,
    int weight_comp)
{
    for (MFIter mfi(input, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto output_arr = output.array(mfi);
        auto input_arr = input.const_array(mfi);
        auto weight_arr = weights.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            output_arr(i,j,k) += weight_arr(i,j,k,weight_comp) * input_arr(i,j,k);
        });
    }
}

inline void accumulateIntoField(MultiFab& output, const MultiFab& input);

inline void accumulateNonlinearVerticalVelocityTerms(
    VerticalVelocityAccumulationContext& context)
{
    auto& workspace = context.workspace;

    for (int eta_order = 0;
         eta_order <= context.perturbation_order - 2;
         ++eta_order)
    {
        const int cumulative_component =
            context.perturbation_order - eta_order - 2;
        const int wavenumber_power = eta_order + 2;

        computeVerticalVelocityTermFromPkComponent(
            workspace.expansion_accum,
            cumulative_component,
            workspace.spectral_work,
            context.fft,
            context.dk_x,
            context.dk_y,
            wavenumber_power);

        accumulateWeightedIntoField(
            context.phi_z,
            workspace.eta_pow,
            workspace.spectral_work,
            eta_order);
    }
}

inline void accumulateLinearVerticalVelocityTerm(
    VerticalVelocityAccumulationContext& context)
{
    auto& workspace = context.workspace;

    computeVerticalVelocityTermFromPkComponent(
        workspace.expansion_accum,
        context.perturbation_order - 1,
        workspace.spectral_work,
        context.fft,
        context.dk_x,
        context.dk_y,
        1);

    accumulateIntoField(context.phi_z, workspace.spectral_work);
}

inline void executeVerticalVelocityAccumulation(
    VerticalVelocityAccumulationContext& context)
{
    auto& workspace = context.workspace;

    initializeCumulativePhiExpansionState(
        workspace.expansion_accum,
        workspace.phi_expansion,
        context.perturbation_order);

    // DNO accumulation must start from zero for both linear and nonlinear cases.
    context.phi_z.setVal(0.0);

    if (context.perturbation_order > 1) {
        accumulateNonlinearVerticalVelocityTerms(context);
    }

    accumulateLinearVerticalVelocityTerm(context);
}

inline void executeVerticalVelocityAccumulation(
    MultiFab& phi_z,
    VerticalVelocityWorkspaceView workspace,
    int perturbation_order,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    Real dk_x,
    Real dk_y)
{
    auto context = makeVerticalVelocityAccumulationContext(
        phi_z,
        workspace,
        perturbation_order,
        fft,
        dk_x,
        dk_y);
    executeVerticalVelocityAccumulation(context);
}

inline void accumulateIntoField(MultiFab& output, const MultiFab& input)
{
    for (MFIter mfi(input, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto output_arr = output.array(mfi);
        auto input_arr = input.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            output_arr(i,j,k) += input_arr(i,j,k);
        });
    }
}

} // namespace crescent::dno::detail
