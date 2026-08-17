// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef FLOODMAR_CHOPPED_NEWTON_SOLVER_HH
#define FLOODMAR_CHOPPED_NEWTON_SOLVER_HH

#include <algorithm>
#include <cstddef>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/assembly/partialreassembler.hh>
#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

/*!
 * \brief Newton solver that keeps 2pnc iterates inside a physically useful
 *        neighborhood before the primary-variable switch is evaluated.
 *
 * The drywell failure is preceded by hypothetical gas-phase mole-fraction
 * sums of O(10)-O(1000).  Those values are Newton overshoots, not a physical
 * phase-boundary crossing.  This solver follows DuMuX's documented chopped
 * update hook and limits pressure, saturation, and composition changes before
 * solutionChanged_ updates volume variables and invokes phase switching.
 */
template<class Assembler, class LinearSolver,
         class Reassembler = PartialReassembler<Assembler>,
         class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator>>
class FloodMarChoppedNewtonSolver
: public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using Scalar = typename Assembler::Scalar;
    using ParentType =
        NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    using VolumeVariables =
        typename Assembler::GridVariables::VolumeVariables;
    using Indices = typename VolumeVariables::Indices;

    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int switchIdx = Indices::switchIdx;
    static constexpr int oxygenIdx = Indices::switchIdx + 1;

    using typename ParentType::Backend;
    using typename ParentType::SolutionVector;
    using typename ParentType::ResidualVector;

public:
    using ParentType::ParentType;
    using typename ParentType::Variables;

private:
    void choppedUpdate_(Variables& varsCurrentIter,
                        const SolutionVector& uLastIter,
                        const ResidualVector& deltaU) final
    {
        auto uCurrentIter = uLastIter;
        Backend::axpy(-1.0, deltaU, uCurrentIter);

        const Scalar maxPressureChange =
            getParam<Scalar>("Newton.ChopMaxPressureChange", 2.0e4);
        const Scalar maxSaturationChange =
            getParam<Scalar>("Newton.ChopMaxSaturationChange", 0.10);
        const Scalar maxLiquidMoleFractionChange =
            getParam<Scalar>("Newton.ChopMaxLiquidMoleFractionChange", 1.0e-5);
        const Scalar maxGasMoleFractionChange =
            getParam<Scalar>("Newton.ChopMaxGasMoleFractionChange", 0.05);
        const Scalar saturationBuffer =
            getParam<Scalar>("Newton.ChopSaturationBuffer", 0.01);

        constexpr Scalar moleFractionFloor = 0.0;
        constexpr Scalar moleFractionCeiling = 1.0 - 1.0e-12;

        const auto boundedStep = [](Scalar candidate,
                                    Scalar oldValue,
                                    Scalar maxChange)
        {
            return std::clamp(
                candidate,
                oldValue - maxChange,
                oldValue + maxChange
            );
        };

        for (std::size_t dofIdx = 0;
             dofIdx < uCurrentIter.size(); ++dofIdx)
        {
            const int phasePresence = uLastIter[dofIdx].state();

            // Avoid pressure jumps of several water-column metres within a
            // single Newton iteration near the moving drywell water surface.
            uCurrentIter[dofIdx][pressureIdx] = boundedStep(
                uCurrentIter[dofIdx][pressureIdx],
                uLastIter[dofIdx][pressureIdx],
                maxPressureChange
            );

            if (phasePresence == Indices::bothPhases)
            {
                // In p0s1, the switching variable is gas saturation.
                uCurrentIter[dofIdx][switchIdx] = std::clamp(
                    boundedStep(
                        uCurrentIter[dofIdx][switchIdx],
                        uLastIter[dofIdx][switchIdx],
                        maxSaturationChange
                    ),
                    -saturationBuffer,
                    Scalar(1.0) + saturationBuffer
                );

                // The additional-component variable is the liquid-phase O2
                // mole fraction and must not jump by orders of magnitude.
                uCurrentIter[dofIdx][oxygenIdx] = std::clamp(
                    boundedStep(
                        uCurrentIter[dofIdx][oxygenIdx],
                        uLastIter[dofIdx][oxygenIdx],
                        maxLiquidMoleFractionChange
                    ),
                    moleFractionFloor,
                    moleFractionCeiling
                );
            }
            else if (phasePresence == Indices::firstPhaseOnly)
            {
                // In liquid-only state, switchIdx and oxygenIdx are dissolved
                // N2 and O2 mole fractions.  Keep each update local and keep
                // their sum below one so the inferred H2O fraction is valid.
                Scalar xN2 = std::clamp(
                    boundedStep(
                        uCurrentIter[dofIdx][switchIdx],
                        uLastIter[dofIdx][switchIdx],
                        maxLiquidMoleFractionChange
                    ),
                    moleFractionFloor,
                    moleFractionCeiling
                );
                Scalar xO2 = std::clamp(
                    boundedStep(
                        uCurrentIter[dofIdx][oxygenIdx],
                        uLastIter[dofIdx][oxygenIdx],
                        maxLiquidMoleFractionChange
                    ),
                    moleFractionFloor,
                    moleFractionCeiling
                );

                const Scalar sum = xN2 + xO2;
                if (sum > moleFractionCeiling)
                {
                    const Scalar scale = moleFractionCeiling/sum;
                    xN2 *= scale;
                    xO2 *= scale;
                }

                uCurrentIter[dofIdx][switchIdx] = xN2;
                uCurrentIter[dofIdx][oxygenIdx] = xO2;
            }
            else if (phasePresence == Indices::secondPhaseOnly)
            {
                // In gas-only state the composition variables can vary over
                // a wider range, but still have to form a valid mixture.
                Scalar xH2O = std::clamp(
                    boundedStep(
                        uCurrentIter[dofIdx][switchIdx],
                        uLastIter[dofIdx][switchIdx],
                        maxGasMoleFractionChange
                    ),
                    moleFractionFloor,
                    moleFractionCeiling
                );
                Scalar xO2 = std::clamp(
                    boundedStep(
                        uCurrentIter[dofIdx][oxygenIdx],
                        uLastIter[dofIdx][oxygenIdx],
                        maxGasMoleFractionChange
                    ),
                    moleFractionFloor,
                    moleFractionCeiling
                );

                const Scalar sum = xH2O + xO2;
                if (sum > moleFractionCeiling)
                {
                    const Scalar scale = moleFractionCeiling/sum;
                    xH2O *= scale;
                    xO2 *= scale;
                }

                uCurrentIter[dofIdx][switchIdx] = xH2O;
                uCurrentIter[dofIdx][oxygenIdx] = xO2;
            }
        }

        this->solutionChanged_(varsCurrentIter, uCurrentIter);

        if (this->enableResidualCriterion())
            this->computeResidualReduction_(varsCurrentIter);
    }
};

} // namespace Dumux

#endif
