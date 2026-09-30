// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef FLOODMAR_PRIMARY_VARIABLE_SWITCH_PHASEGUARD_V2_HH
#define FLOODMAR_PRIMARY_VARIABLE_SWITCH_PHASEGUARD_V2_HH

#include <algorithm>
#include <iostream>
#include <vector>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \brief 2pnc primary-variable switch for FloodMAR.
 *
 * Design goals:
 *
 * 1. Preserve a true thermodynamic appearance hysteresis:
 *       gas appears only if hypothetical gas mole-fraction sum
 *       exceeds 1 + GasAppearanceTolerance.
 *
 * 2. Remove the iteration-count gas hold used by RestartV4.
 *
 * 3. If a newly appearing gas phase becomes non-physical
 *       (Sgas <= 0)
 *    during the same nonlinear solve, immediately remove it.
 *
 * 4. Once gas disappears at a DOF, prohibit immediate
 *    reappearance during the same nonlinear solve.
 *
 * reset() clears the per-solve latch. Therefore a later
 * nonlinear solve / reduced-dt retry is free to determine
 * the physically appropriate phase state again.
 */
class FloodMarPhaseGuardV2PrimaryVariableSwitch
: public PrimaryVariableSwitch<FloodMarPhaseGuardV2PrimaryVariableSwitch>
{
    using ParentType =
        PrimaryVariableSwitch<FloodMarPhaseGuardV2PrimaryVariableSwitch>;

    friend ParentType;

public:

    explicit FloodMarPhaseGuardV2PrimaryVariableSwitch(int verbosity = 1)
    : ParentType(verbosity)
    , gasAppearanceTolerance_(
        std::clamp(
            getParam<double>(
                "PrimaryVariableSwitch.GasAppearanceTolerance",
                0.02
            ),
            0.0,
            0.20
        )
      )
    , liquidAppearanceTolerance_(
        std::clamp(
            getParam<double>(
                "PrimaryVariableSwitch.LiquidAppearanceTolerance",
                0.02
            ),
            0.0,
            0.20
        )
      )
    , appearanceSaturation_(
        std::clamp(
            getParam<double>(
                "PrimaryVariableSwitch.AppearanceSaturation",
                1.0e-4
            ),
            1.0e-12,
            1.0e-2
        )
      )
    {}

    void reset(const std::size_t numDofs)
    {
        ParentType::reset(numDofs);

        //  0 : no gas transition yet in this nonlinear solve
        // +1 : gas appeared
        // -1 : gas disappeared; block reappearance until reset()
        gasTransitionState_.assign(numDofs, 0);
    }

protected:

    template<class VolumeVariables,
             class IndexType,
             class GlobalPosition>
    bool update_(
        typename VolumeVariables::PrimaryVariables& priVars,
        const VolumeVariables& volVars,
        IndexType dofIdxGlobal,
        const GlobalPosition& globalPos
    )
    {
        using Scalar =
            typename VolumeVariables::PrimaryVariables::value_type;

        using FluidSystem =
            typename VolumeVariables::FluidSystem;

        using Indices =
            typename VolumeVariables::Indices;

        static constexpr int phase0Idx =
            FluidSystem::phase0Idx;

        static constexpr int phase1Idx =
            FluidSystem::phase1Idx;

        static constexpr int comp0Idx =
            FluidSystem::comp0Idx;

        static constexpr int comp1Idx =
            FluidSystem::comp1Idx;

        static constexpr auto numComponents =
            VolumeVariables::numFluidComponents();

        static constexpr auto numMajorComponents =
            VolumeVariables::numFluidPhases();

        static constexpr bool useMoles =
            VolumeVariables::useMoles();

        static constexpr auto formulation =
            VolumeVariables::priVarFormulation();

        static constexpr int switchIdx =
            Indices::switchIdx;

        static_assert(
            formulation == TwoPFormulation::p0s1
            || formulation == TwoPFormulation::p1s0,
            "Chosen TwoPFormulation is not supported"
        );

        static_assert(
            useMoles || numComponents < 3,
            "Mass-fraction formulation is only implemented "
            "for fewer than three components"
        );

        const std::size_t dofIdx =
            static_cast<std::size_t>(dofIdxGlobal);

        const int phasePresence =
            priVars.state();

        int newPhasePresence =
            phasePresence;

        bool switched = false;


        // ==============================================================
        // BOTH PHASES PRESENT
        // ==============================================================

        if (phasePresence == Indices::bothPhases)
        {
            // ----------------------------------------------------------
            // Phase 0 disappearance:
            // retain the established DuMuX/RestartV4 behavior.
            // ----------------------------------------------------------

            Scalar phase0Smin = 0.0;

            if (this->wasSwitched_[dofIdx])
                phase0Smin = -0.01;

            if (volVars.saturation(phase0Idx) <= phase0Smin)
            {
                switched = true;

                if (this->verbosity() > 1)
                    std::cout
                        << "First phase ("
                        << FluidSystem::phaseName(phase0Idx)
                        << ") disappears at dof "
                        << dofIdxGlobal
                        << ", coordinates: "
                        << globalPos
                        << ", S_"
                        << FluidSystem::phaseName(phase0Idx)
                        << ": "
                        << volVars.saturation(phase0Idx)
                        << std::endl;

                newPhasePresence =
                    Indices::secondPhaseOnly;

                if constexpr (useMoles)
                    priVars[switchIdx] =
                        volVars.moleFraction(
                            phase1Idx,
                            comp0Idx
                        );
                else
                    priVars[switchIdx] =
                        volVars.massFraction(
                            phase1Idx,
                            comp0Idx
                        );

                if constexpr (useMoles)
                {
                    for (int compIdx = numMajorComponents;
                         compIdx < numComponents;
                         ++compIdx)
                    {
                        priVars[compIdx] =
                            volVars.moleFraction(
                                phase1Idx,
                                compIdx
                            );
                    }
                }
                else
                {
                    for (int compIdx = numMajorComponents;
                         compIdx < numComponents;
                         ++compIdx)
                    {
                        priVars[compIdx] =
                            volVars.massFraction(
                                phase1Idx,
                                compIdx
                            );
                    }
                }
            }

            // ----------------------------------------------------------
            // Gas-phase disappearance.
            //
            // KEY CHANGE:
            //
            // No 48-iteration hold.
            // No accepted degenerate zero-gas two-phase state.
            //
            // As soon as Sgas <= 0, return to liquid-only and latch
            // gas off for the remainder of this nonlinear solve.
            // ----------------------------------------------------------

            else if (volVars.saturation(phase1Idx) <= Scalar(0.0))
            {
                switched = true;

                gasTransitionState_[dofIdx] = -1;

                if (this->verbosity() > 1)
                    std::cout
                        << "Second phase ("
                        << FluidSystem::phaseName(phase1Idx)
                        << ") disappears at dof "
                        << dofIdxGlobal
                        << ", coordinates: "
                        << globalPos
                        << ", S_"
                        << FluidSystem::phaseName(phase1Idx)
                        << ": "
                        << volVars.saturation(phase1Idx)
                        << " [PhaseGuardV2 latch]"
                        << std::endl;

                newPhasePresence =
                    Indices::firstPhaseOnly;

                if constexpr (useMoles)
                    priVars[switchIdx] =
                        volVars.moleFraction(
                            phase0Idx,
                            comp1Idx
                        );
                else
                    priVars[switchIdx] =
                        volVars.massFraction(
                            phase0Idx,
                            comp1Idx
                        );
            }
        }


        // ==============================================================
        // PHASE 1 ONLY
        // ==============================================================

        else if (phasePresence == Indices::secondPhaseOnly)
        {
            Scalar x0Sum = 0.0;

            for (int compIdx = 0;
                 compIdx < numComponents;
                 ++compIdx)
            {
                x0Sum +=
                    volVars.moleFraction(
                        phase0Idx,
                        compIdx
                    );
            }

            const Scalar threshold =
                Scalar(1.0)
                + static_cast<Scalar>(
                    liquidAppearanceTolerance_
                  );

            if (x0Sum > threshold)
            {
                switched = true;

                if (this->verbosity() > 1)
                    std::cout
                        << "First phase ("
                        << FluidSystem::phaseName(phase0Idx)
                        << ") appears at dof "
                        << dofIdxGlobal
                        << ", coordinates: "
                        << globalPos
                        << ", sum x^i_"
                        << FluidSystem::phaseName(phase0Idx)
                        << ": "
                        << x0Sum
                        << ", threshold: "
                        << threshold
                        << std::endl;

                newPhasePresence =
                    Indices::bothPhases;

                priVars[switchIdx] =
                    formulation == TwoPFormulation::p1s0
                    ? static_cast<Scalar>(
                        appearanceSaturation_
                      )
                    : Scalar(1.0)
                      - static_cast<Scalar>(
                          appearanceSaturation_
                        );

                for (int compIdx = numMajorComponents;
                     compIdx < numComponents;
                     ++compIdx)
                {
                    priVars[compIdx] =
                        volVars.moleFraction(
                            phase0Idx,
                            compIdx
                        );
                }
            }
        }


        // ==============================================================
        // PHASE 0 / LIQUID ONLY
        // ==============================================================

        else if (phasePresence == Indices::firstPhaseOnly)
        {
            Scalar x1Sum = 0.0;

            for (int compIdx = 0;
                 compIdx < numComponents;
                 ++compIdx)
            {
                x1Sum +=
                    volVars.moleFraction(
                        phase1Idx,
                        compIdx
                    );
            }

            const Scalar threshold =
                Scalar(1.0)
                + static_cast<Scalar>(
                    gasAppearanceTolerance_
                  );

            // If gas disappeared earlier in THIS nonlinear solve,
            // do not let the same DOF immediately create it again.
            if (gasTransitionState_[dofIdx] >= 0
                && x1Sum > threshold)
            {
                switched = true;

                if (this->verbosity() > 1)
                    std::cout
                        << "Second phase ("
                        << FluidSystem::phaseName(phase1Idx)
                        << ") appears at dof "
                        << dofIdxGlobal
                        << ", coordinates: "
                        << globalPos
                        << ", sum x^i_"
                        << FluidSystem::phaseName(phase1Idx)
                        << ": "
                        << x1Sum
                        << ", threshold: "
                        << threshold
                        << " [PhaseGuardV2]"
                        << std::endl;

                newPhasePresence =
                    Indices::bothPhases;

                priVars[switchIdx] =
                    formulation == TwoPFormulation::p1s0
                    ? Scalar(1.0)
                      - static_cast<Scalar>(
                          appearanceSaturation_
                        )
                    : static_cast<Scalar>(
                        appearanceSaturation_
                      );

                gasTransitionState_[dofIdx] = 1;
            }
        }


        priVars.setState(newPhasePresence);

        this->wasSwitched_[dofIdx] =
            switched;

        return phasePresence !=
               newPhasePresence;
    }


private:

    double gasAppearanceTolerance_;
    double liquidAppearanceTolerance_;
    double appearanceSaturation_;

    std::vector<signed char>
        gasTransitionState_;
};


template<class Traits>
class FloodMarPhaseGuardV2TwoPNCVolumeVariables
: public TwoPNCVolumeVariables<Traits>
{
    using ParentType =
        TwoPNCVolumeVariables<Traits>;

public:

    using ParentType::ParentType;

    using PrimaryVariableSwitch =
        FloodMarPhaseGuardV2PrimaryVariableSwitch;
};

} // namespace Dumux

#endif
