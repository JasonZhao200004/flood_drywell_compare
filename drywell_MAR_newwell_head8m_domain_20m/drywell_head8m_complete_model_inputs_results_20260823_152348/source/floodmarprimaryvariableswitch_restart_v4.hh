// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef FLOODMAR_PRIMARY_VARIABLE_SWITCH_RESTART_V4_HH
#define FLOODMAR_PRIMARY_VARIABLE_SWITCH_RESTART_V4_HH

#include <algorithm>
#include <iostream>
#include <vector>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \brief 2pnc primary-variable switch with a true appearance/disappearance
 *        hysteresis gap and an intra-solve gas hold.
 *
 * A liquid-only cell does not create gas until the hypothetical gas-phase
 * mole-fraction sum exceeds 1 + GasAppearanceTolerance.  A two-phase cell
 * still loses gas at zero saturation unless gas appeared during the current
 * nonlinear solve, in which case the existing negative-saturation hold is
 * used.  The gap prevents the same accepted state from triggering an
 * identical active-set jump after every time-step cut.
 */
class FloodMarRestartV4PrimaryVariableSwitch
: public PrimaryVariableSwitch<FloodMarRestartV4PrimaryVariableSwitch>
{
    using ParentType =
        PrimaryVariableSwitch<FloodMarRestartV4PrimaryVariableSwitch>;

    friend ParentType;

public:
    explicit FloodMarRestartV4PrimaryVariableSwitch(int verbosity = 1)
    : ParentType(verbosity)
    , gasPhaseHoldIterations_(
        std::max(1, getParam<int>(
            "PrimaryVariableSwitch.GasPhaseHoldIterations", 48)))
    , gasDisappearanceTolerance_(
        std::max(0.0, getParam<double>(
            "PrimaryVariableSwitch.GasDisappearanceTolerance", 0.01)))
    , gasAppearanceTolerance_(
        std::clamp(getParam<double>(
            "PrimaryVariableSwitch.GasAppearanceTolerance", 0.02),
            0.0,
            0.20))
    , liquidAppearanceTolerance_(
        std::clamp(getParam<double>(
            "PrimaryVariableSwitch.LiquidAppearanceTolerance", 0.02),
            0.0,
            0.20))
    , appearanceSaturation_(
        std::clamp(getParam<double>(
            "PrimaryVariableSwitch.AppearanceSaturation", 1.0e-4),
            1.0e-12,
            1.0e-2))
    {}

    void reset(const std::size_t numDofs)
    {
        ParentType::reset(numDofs);
        gasPhaseHoldRemaining_.assign(numDofs, 0);
    }

protected:
    template<class VolumeVariables, class IndexType, class GlobalPosition>
    bool update_(typename VolumeVariables::PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 IndexType dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        using Scalar =
            typename VolumeVariables::PrimaryVariables::value_type;
        using FluidSystem = typename VolumeVariables::FluidSystem;
        using Indices = typename VolumeVariables::Indices;

        static constexpr int phase0Idx = FluidSystem::phase0Idx;
        static constexpr int phase1Idx = FluidSystem::phase1Idx;
        static constexpr int comp0Idx = FluidSystem::comp0Idx;
        static constexpr int comp1Idx = FluidSystem::comp1Idx;
        static constexpr auto numComponents =
            VolumeVariables::numFluidComponents();
        static constexpr auto numMajorComponents =
            VolumeVariables::numFluidPhases();
        static constexpr bool useMoles = VolumeVariables::useMoles();
        static constexpr auto formulation =
            VolumeVariables::priVarFormulation();
        static constexpr int switchIdx = Indices::switchIdx;

        static_assert(
            formulation == TwoPFormulation::p0s1
            || formulation == TwoPFormulation::p1s0,
            "Chosen TwoPFormulation is not supported");
        static_assert(
            useMoles || numComponents < 3,
            "Mass-fraction formulation is only implemented for fewer than three components");

        const std::size_t dofIdx =
            static_cast<std::size_t>(dofIdxGlobal);

        const int phasePresence = priVars.state();
        int newPhasePresence = phasePresence;
        bool switched = false;

        if (phasePresence == Indices::bothPhases)
        {
            Scalar phase0Smin = 0.0;
            if (this->wasSwitched_[dofIdx])
                phase0Smin = -0.01;

            if (volVars.saturation(phase0Idx) <= phase0Smin)
            {
                switched = true;
                if (this->verbosity() > 1)
                    std::cout
                        << "First phase (" << FluidSystem::phaseName(phase0Idx)
                        << ") disappears at dof " << dofIdxGlobal
                        << ", coordinates: " << globalPos
                        << ", S_" << FluidSystem::phaseName(phase0Idx)
                        << ": " << volVars.saturation(phase0Idx)
                        << std::endl;

                newPhasePresence = Indices::secondPhaseOnly;
                if constexpr (useMoles)
                    priVars[switchIdx] =
                        volVars.moleFraction(phase1Idx, comp0Idx);
                else
                    priVars[switchIdx] =
                        volVars.massFraction(phase1Idx, comp0Idx);

                if constexpr (useMoles)
                    for (int compIdx = numMajorComponents;
                         compIdx < numComponents; ++compIdx)
                        priVars[compIdx] =
                            volVars.moleFraction(phase1Idx, compIdx);
                else
                    for (int compIdx = numMajorComponents;
                         compIdx < numComponents; ++compIdx)
                        priVars[compIdx] =
                            volVars.massFraction(phase1Idx, compIdx);
            }
            else
            {
                Scalar gasSmin = 0.0;
                if (gasPhaseHoldRemaining_[dofIdx] > 0)
                {
                    gasSmin =
                        -static_cast<Scalar>(gasDisappearanceTolerance_);
                    --gasPhaseHoldRemaining_[dofIdx];
                }
                else if (this->wasSwitched_[dofIdx])
                    gasSmin = -0.01;

                if (volVars.saturation(phase1Idx) <= gasSmin)
                {
                    switched = true;
                    gasPhaseHoldRemaining_[dofIdx] = 0;

                    if (this->verbosity() > 1)
                        std::cout
                            << "Second phase (" << FluidSystem::phaseName(phase1Idx)
                            << ") disappears at dof " << dofIdxGlobal
                            << ", coordinates: " << globalPos
                            << ", S_" << FluidSystem::phaseName(phase1Idx)
                            << ": " << volVars.saturation(phase1Idx)
                            << std::endl;

                    newPhasePresence = Indices::firstPhaseOnly;
                    if constexpr (useMoles)
                        priVars[switchIdx] =
                            volVars.moleFraction(phase0Idx, comp1Idx);
                    else
                        priVars[switchIdx] =
                            volVars.massFraction(phase0Idx, comp1Idx);
                }
            }
        }
        else if (phasePresence == Indices::secondPhaseOnly)
        {
            Scalar x0Sum = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                x0Sum += volVars.moleFraction(phase0Idx, compIdx);

            const Scalar threshold =
                Scalar(1.0)
                + static_cast<Scalar>(liquidAppearanceTolerance_);

            if (x0Sum > threshold)
            {
                switched = true;
                if (this->verbosity() > 1)
                    std::cout
                        << "First phase (" << FluidSystem::phaseName(phase0Idx)
                        << ") appears at dof " << dofIdxGlobal
                        << ", coordinates: " << globalPos
                        << ", sum x^i_" << FluidSystem::phaseName(phase0Idx)
                        << ": " << x0Sum
                        << ", threshold: " << threshold
                        << std::endl;

                newPhasePresence = Indices::bothPhases;
                priVars[switchIdx] =
                    formulation == TwoPFormulation::p1s0
                    ? static_cast<Scalar>(appearanceSaturation_)
                    : Scalar(1.0) - static_cast<Scalar>(appearanceSaturation_);

                for (int compIdx = numMajorComponents;
                     compIdx < numComponents; ++compIdx)
                    priVars[compIdx] =
                        volVars.moleFraction(phase0Idx, compIdx);
            }
        }
        else if (phasePresence == Indices::firstPhaseOnly)
        {
            Scalar x1Sum = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                x1Sum += volVars.moleFraction(phase1Idx, compIdx);

            const Scalar threshold =
                Scalar(1.0)
                + static_cast<Scalar>(gasAppearanceTolerance_);

            if (x1Sum > threshold)
            {
                switched = true;
                if (this->verbosity() > 1)
                    std::cout
                        << "Second phase (" << FluidSystem::phaseName(phase1Idx)
                        << ") appears at dof " << dofIdxGlobal
                        << ", coordinates: " << globalPos
                        << ", sum x^i_" << FluidSystem::phaseName(phase1Idx)
                        << ": " << x1Sum
                        << ", threshold: " << threshold
                        << std::endl;

                newPhasePresence = Indices::bothPhases;
                priVars[switchIdx] =
                    formulation == TwoPFormulation::p1s0
                    ? Scalar(1.0) - static_cast<Scalar>(appearanceSaturation_)
                    : static_cast<Scalar>(appearanceSaturation_);
                gasPhaseHoldRemaining_[dofIdx] =
                    gasPhaseHoldIterations_;
            }
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdx] = switched;
        return phasePresence != newPhasePresence;
    }

private:
    int gasPhaseHoldIterations_;
    double gasDisappearanceTolerance_;
    double gasAppearanceTolerance_;
    double liquidAppearanceTolerance_;
    double appearanceSaturation_;
    std::vector<int> gasPhaseHoldRemaining_;
};

template<class Traits>
class FloodMarRestartV4TwoPNCVolumeVariables
: public TwoPNCVolumeVariables<Traits>
{
    using ParentType = TwoPNCVolumeVariables<Traits>;

public:
    using ParentType::ParentType;
    using PrimaryVariableSwitch =
        FloodMarRestartV4PrimaryVariableSwitch;
};

} // namespace Dumux

#endif
