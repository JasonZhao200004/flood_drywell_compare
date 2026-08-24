// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// SPDX-License-Identifier: GPL-3.0-or-later

#include <config.h>

#include <algorithm>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/assembly/fvassembler.hh>

// This model uses an unstructured UGGrid loaded from floodmar.msh.
#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;

    // Initialize MPI and the selected multithreading backend.
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // Read the input file and command-line overrides.
    Parameters::init(argc, argv);

    // Read the unstructured Gmsh mesh with UGGrid.
    using Grid = GetPropType<TypeTag, Properties::Grid>;

    GridManager<Grid> gridManager;
    gridManager.init();

    const auto& leafGridView =
        gridManager.grid().leafGridView();

    if (mpiHelper.rank() == 0)
    {
        std::cout
            << "UGGrid loaded successfully: "
            << leafGridView.size(0)
            << " elements and "
            << leafGridView.size(Grid::dimension)
            << " vertices."
            << std::endl;
    }

    // Construct the axisymmetric finite-volume grid geometry.
    using GridGeometry =
        GetPropType<TypeTag, Properties::GridGeometry>;

    auto gridGeometry =
        std::make_shared<GridGeometry>(
            leafGridView
        );

    // Create the problem containing initial and boundary conditions.
    using Problem =
        GetPropType<TypeTag, Properties::Problem>;

    auto problem =
        std::make_shared<Problem>(
            gridGeometry
        );

    // Initialize the solution vector.
    using SolutionVector =
        GetPropType<TypeTag, Properties::SolutionVector>;

    SolutionVector x;
    problem->applyInitialSolution(x);

    auto xOld = x;

    // Initialize volume and flux variables.
    using GridVariables =
        GetPropType<TypeTag, Properties::GridVariables>;

    auto gridVariables =
        std::make_shared<GridVariables>(
            problem,
            gridGeometry
        );

    gridVariables->init(x);

    // Configure VTK output for ParaView.
    VtkOutputModule<GridVariables, SolutionVector>
        vtkWriter(
            *gridVariables,
            x,
            problem->name()
        );

    using VelocityOutput =
        GetPropType<TypeTag, Properties::VelocityOutput>;

    vtkWriter.addVelocityOutput(
        std::make_shared<VelocityOutput>(
            *gridVariables
        )
    );

    GetPropType<
        TypeTag,
        Properties::IOFields
    >::initOutputModule(vtkWriter);

    // Write the initial condition.
    vtkWriter.write(0.0);

    // Read time-loop parameters. DuMuX uses seconds here.
    using Scalar =
        GetPropType<TypeTag, Properties::Scalar>;

    const Scalar tEnd =
        getParam<Scalar>(
            "TimeLoop.TEnd"
        );

    const Scalar dtInitial =
        getParam<Scalar>(
            "TimeLoop.DtInitial"
        );

    const Scalar maxDt =
        getParam<Scalar>(
            "TimeLoop.MaxTimeStepSize"
        );

    const Scalar outputInterval =
        getParam<Scalar>(
            "TimeLoop.OutputInterval",
            21600.0
        );

    auto timeLoop =
        std::make_shared<TimeLoop<Scalar>>(
            0.0,
            dtInitial,
            tEnd
        );

    timeLoop->setMaxTimeStepSize(maxDt);

    // Create the finite-volume assembler.
    using Assembler =
        FVAssembler<TypeTag, DiffMethod::numeric>;

    auto assembler =
        std::make_shared<Assembler>(
            problem,
            gridGeometry,
            gridVariables,
            timeLoop,
            xOld
        );

    // Create the linear solver.
    using LinearSolver =
        AMGBiCGSTABIstlSolver<
            LinearSolverTraits<GridGeometry>,
            LinearAlgebraTraitsFromAssembler<Assembler>
        >;

    auto linearSolver =
        std::make_shared<LinearSolver>(
            leafGridView,
            gridGeometry->dofMapper()
        );

    // Create the nonlinear Newton solver.
    using NonlinearSolver =
        NewtonSolver<Assembler, LinearSolver>;

    NonlinearSolver nonlinearSolver(
        assembler,
        linearSolver
    );

    Scalar nextOutputTime =
        outputInterval;

    timeLoop->start();

    do
    {
        // Never step across a change point in the time-variable drywell
        // boundary. This makes the 5-hour linear ramps reproducible.
        const Scalar nextBoundaryEvent =
            problem->nextBoundaryEventTime(timeLoop->time());
        if (nextBoundaryEvent < timeLoop->time() + timeLoop->timeStepSize())
        {
            timeLoop->setTimeStepSize(
                nextBoundaryEvent - timeLoop->time()
            );
        }

        // Evaluate the time-variable drywell boundary condition
        // at the new implicit time level.
        problem->setTime(
            timeLoop->time()
            + timeLoop->timeStepSize()
        );

        // Solve the nonlinear two-phase system.
        nonlinearSolver.solve(
            x,
            *timeLoop
        );

        // Accept the new solution.
        xOld = x;
        gridVariables->advanceTimeStep();

        // Advance simulation time.
        timeLoop->advanceTimeStep();
        timeLoop->reportTimeStep();

        // Write output at the specified interval and final time.
        const Scalar timeTolerance =
            1.0e-8
            * std::max(
                Scalar(1.0),
                timeLoop->time()
            );

        if (
            timeLoop->time()
                >= nextOutputTime - timeTolerance
            || timeLoop->finished()
        )
        {
            vtkWriter.write(
                timeLoop->time()
            );

            while (
                nextOutputTime
                <= timeLoop->time()
                    + timeTolerance
            )
            {
                nextOutputTime +=
                    outputInterval;
            }
        }

        // Use the Newton solver's suggested next time step,
        // but do not exceed the configured maximum.
        const Scalar suggestedDt =
            nonlinearSolver.suggestTimeStepSize(
                timeLoop->timeStepSize()
            );

        timeLoop->setTimeStepSize(
            std::min(
                maxDt,
                suggestedDt
            )
        );

    } while (!timeLoop->finished());

    timeLoop->finalize(
        leafGridView.comm()
    );

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
