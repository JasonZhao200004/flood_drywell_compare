// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef FLOODMAR_PRIMARY_VARIABLE_SWITCH_PHASEGUARD_HH
#define FLOODMAR_PRIMARY_VARIABLE_SWITCH_PHASEGUARD_HH

#include <algorithm>
#include <iostream>
#include <vector>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \brief 2pnc primary-variable switch with a per-Newton-solve gas-phase guard.
 *
 * The standard active-set switch may change one degree of freedom repeatedly
 * between liquid-only and two-phase states during a difficult Newton solve.
 * This variant permits gas appearance normally, protects a newly appearing
 * gas phase down to a configurable small negative saturation, and then blocks
 * immediate reappearance if it has to disappear again in the same nonlinear
 * solve.  reset() clears the guard before the next nonlinear solve/retry, so
 * the converged phase criterion is reconsidered after the time step is cut.
 */
class FloodMarPhaseGuardPrimaryVariableSwitch
: public PrimaryVariableSwitch<FloodMarPhaseGuardPrimaryVariableSwitch>
{
    using ParentType =
        PrimaryVariableSwitch<FloodMarPhaseGuardPrimaryVariableSwitch>;

    friend ParentType;

public:
    explicit FloodMarPhaseGuardPrimaryVariableSwitch(int verbosity = 1)
    : ParentType(verbosity)
    , gasDisappearanceTolerance_(
        std::max(0.0, getParam<double>(
            "PrimaryVariableSwitch.GasDisappearanceTolerance", 0.01)))
    , gasAppearanceSaturation_(
        std::clamp(getParam<double>(
            "PrimaryVariableSwitch.GasAppearanceSaturation", 1.0e-4),
            1.0e-12,
            1.0e-2))
    {}

    void reset(const std::size_t numDofs)
    {
        ParentType::reset(numDofs);
        gasTransitionState_.assign(numDofs, 0);
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

        bool wouldSwitch = false;
        const int phasePresence = priVars.state();
        int newPhasePresence = phasePresence;

        if (phasePresence == Indices::bothPhases)
        {
            Scalar phase0Smin = 0.0;
            if (this->wasSwitched_[dofIdx])
                phase0Smin = -0.01;

            // Keep the standard rule for disappearance of phase 0.
            if (volVars.saturation(phase0Idx) <= phase0Smin)
            {
                wouldSwitch = true;
                if (this->verbosity() > 1)
                    std::cout
                        << "First phase ("
                        << FluidSystem::phaseName(phase0Idx)
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
                // State +1 means that gas appeared during this nonlinear
                // solve.  Protect it from an immediate numerical reversal,
                // but retain a hard negative-saturation escape threshold.
                Scalar gasSmin = 0.0;
                if (gasTransitionState_[dofIdx] > 0)
                    gasSmin =
                        -static_cast<Scalar>(gasDisappearanceTolerance_);
                else if (this->wasSwitched_[dofIdx])
                    gasSmin = -0.01;

                if (volVars.saturation(phase1Idx) <= gasSmin)
                {
                    wouldSwitch = true;
                    // Once gas disappears, do not let the same DOF bring it
                    // back during this nonlinear solve.  A smaller-dt retry
                    // calls reset() and is free to choose the phase again.
                    gasTransitionState_[dofIdx] = -1;

                    if (this->verbosity() > 1)
                        std::cout
                            << "Second phase ("
                            << FluidSystem::phaseName(phase1Idx)
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
            Scalar x0Max = 1.0;
            Scalar x0Sum = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                x0Sum += volVars.moleFraction(phase0Idx, compIdx);

            if (x0Sum > 1.0)
            {
                wouldSwitch = true;
                if (this->wasSwitched_[dofIdx])
                    x0Max *= 1.02;

                if (x0Sum > x0Max)
                {
                    if (this->verbosity() > 1)
                        std::cout
                            << "First phase ("
                            << FluidSystem::phaseName(phase0Idx)
                            << ") appears at dof " << dofIdxGlobal
                            << ", coordinates: " << globalPos
                            << ", sum x^i_"
                            << FluidSystem::phaseName(phase0Idx)
                            << ": " << x0Sum << std::endl;

                    newPhasePresence = Indices::bothPhases;
                    priVars[switchIdx] =
                        formulation == TwoPFormulation::p1s0
                        ? static_cast<Scalar>(gasAppearanceSaturation_)
                        : Scalar(1.0)
                          - static_cast<Scalar>(gasAppearanceSaturation_);

                    for (int compIdx = numMajorComponents;
                         compIdx < numComponents; ++compIdx)
                        priVars[compIdx] =
                            volVars.moleFraction(phase0Idx, compIdx);
                }
            }
        }
        else if (phasePresence == Indices::firstPhaseOnly)
        {
            Scalar x1Max = 1.0;
            Scalar x1Sum = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                x1Sum += volVars.moleFraction(phase1Idx, compIdx);

            // A gas phase that disappeared in this nonlinear solve is
            // latched off until reset().  This is the key anti-chatter guard.
            if (gasTransitionState_[dofIdx] >= 0 && x1Sum > 1.0)
            {
                wouldSwitch = true;
                if (this->wasSwitched_[dofIdx])
                    x1Max *= 1.02;

                if (x1Sum > x1Max)
                {
                    if (this->verbosity() > 1)
                        std::cout
                            << "Second phase ("
                            << FluidSystem::phaseName(phase1Idx)
                            << ") appears at dof " << dofIdxGlobal
                            << ", coordinates: " << globalPos
                            << ", sum x^i_"
                            << FluidSystem::phaseName(phase1Idx)
                            << ": " << x1Sum << std::endl;

                    newPhasePresence = Indices::bothPhases;
                    priVars[switchIdx] =
                        formulation == TwoPFormulation::p1s0
                        ? Scalar(1.0)
                          - static_cast<Scalar>(gasAppearanceSaturation_)
                        : static_cast<Scalar>(gasAppearanceSaturation_);
                    gasTransitionState_[dofIdx] = 1;
                }
            }
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdx] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

private:
    double gasDisappearanceTolerance_;
    double gasAppearanceSaturation_;

    // 0: no gas transition in this solve; +1: appeared; -1: disappeared.
    std::vector<signed char> gasTransitionState_;
};

template<class Traits>
class FloodMarPhaseGuardTwoPNCVolumeVariables
: public TwoPNCVolumeVariables<Traits>
{
    using ParentType = TwoPNCVolumeVariables<Traits>;

public:
    using ParentType::ParentType;
    using PrimaryVariableSwitch =
        FloodMarPhaseGuardPrimaryVariableSwitch;
};

} // namespace Dumux

#endif
