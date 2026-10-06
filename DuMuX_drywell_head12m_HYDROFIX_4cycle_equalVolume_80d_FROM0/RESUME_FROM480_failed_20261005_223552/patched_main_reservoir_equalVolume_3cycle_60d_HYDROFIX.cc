// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// SPDX-License-Identifier: GPL-3.0-or-later

#include <config.h>

#include <algorithm>
#include <map>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

#include "properties_reservoir_equalVolume_3cycle_60d_HYDROFIX.hh"
#include "floodmarchoppednewtonsolver_restart_v4.hh"

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

    using Scalar =
        GetPropType<TypeTag, Properties::Scalar>;

    // Initialize the solution vector, optionally replacing the analytical
    // initial condition with a cell-centred FloodMAR restart state.
    using SolutionVector =
        GetPropType<TypeTag, Properties::SolutionVector>;

    SolutionVector x;
    problem->applyInitialSolution(x);

    Scalar startTime = 0.0;
    const bool restartEnabled =
        getParam<bool>("Restart.Enable", false);

    if (restartEnabled)
    {
        const std::string restartFile =
            getParam<std::string>("Restart.File");
        std::ifstream input(restartFile);
        if (!input)
            throw std::runtime_error(
                "Could not open FloodMAR restart file: " + restartFile);

        std::string magic;
        std::size_t numDofs = 0;
        input >> magic >> numDofs >> startTime;

        if (magic != "FLOODMAR_RESTART_V1")
            throw std::runtime_error(
                "Unsupported FloodMAR restart format in: " + restartFile);
        if (numDofs != x.size())
            throw std::runtime_error(
                "Restart DOF count does not match the current grid");

        std::vector<bool> seen(numDofs, false);
        for (std::size_t row = 0; row < numDofs; ++row)
        {
            std::size_t dofIdx = 0;
            int phaseState = 0;
            Scalar pressure = 0.0;
            Scalar switchValue = 0.0;
            Scalar oxygen = 0.0;

            if (!(input
                  >> dofIdx
                  >> phaseState
                  >> pressure
                  >> switchValue
                  >> oxygen))
                throw std::runtime_error(
                    "Restart file ended before all DOFs were read");

            if (dofIdx >= numDofs || seen[dofIdx])
                throw std::runtime_error(
                    "Restart file contains an invalid or duplicate DOF index");

            x[dofIdx].setState(phaseState);
            x[dofIdx][0] = pressure;
            x[dofIdx][1] = switchValue;
            x[dofIdx][2] = oxygen;
            seen[dofIdx] = true;
        }

        std::size_t extraDof = 0;
        if (input >> extraDof)
            throw std::runtime_error(
                "Restart file contains more rows than the current grid");

        problem->setTime(startTime);

        if (mpiHelper.rank() == 0)
            std::cout
                << "Loaded FloodMAR restart: " << restartFile
                << "\nRestart time = " << startTime
                << " s (" << startTime/3600.0 << " h)"
                << std::endl;
    }


    struct RecoveryRow { double cx,cy,phase,pl,sg,n2,o2,pg,og,ng; };
    std::vector<RecoveryRow> recovered(x.size());
    const bool recoveryMode = getParam<bool>("Recovery.Enable", false);
    if (recoveryMode)
    {
        std::ifstream in(getParam<std::string>("Recovery.File"));
        std::size_t n=0;
        if (!(in >> n >> startTime) || n!=x.size())
            throw std::runtime_error("Invalid recovery file DOF count");
        std::map<std::pair<long long,long long>, RecoveryRow> lookup;
        auto key=[](double a,double b) { return std::make_pair(std::llround(a*1e8),std::llround(b*1e8)); };
        for (std::size_t i=0;i<n;++i) {
            RecoveryRow r;
            if (!(in>>r.cx>>r.cy>>r.phase>>r.pl>>r.sg>>r.n2>>r.o2>>r.pg>>r.og>>r.ng))
                throw std::runtime_error("Incomplete recovery cell data");
            if (!lookup.emplace(key(r.cx,r.cy),r).second)
                throw std::runtime_error("Duplicate cell centers in VTU recovery");
        }
        using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
        for (const auto& element : elements(leafGridView)) {
            const auto pos=element.geometry().center();
            auto it=lookup.find(key(pos[0],pos[1]));
            if(it==lookup.end()) throw std::runtime_error("VTU cell center does not match mesh; recovery stopped");
            const auto r=it->second;
            const auto idx=gridGeometry->elementMapper().index(element);
            const int phase=static_cast<int>(std::llround(r.phase));
            if (std::abs(r.phase-phase)>1e-8) throw std::runtime_error("Invalid phase state");
            x[idx].setState(phase);
            x[idx][0]=r.pl;
            if(phase==Indices::bothPhases) x[idx][1]=r.sg;
            else if(phase==Indices::firstPhaseOnly) x[idx][1]=r.n2;
            else throw std::runtime_error("Gas-only cell encountered; cannot safely recover with liquid-pressure formulation");
            x[idx][2]=r.o2;
            recovered[idx]=r;
        }
        problem->setTime(startTime);
        std::cout << "Recovered cell-centered VTU at " << startTime/3600.0 << " h\n";
    }
    if (restartEnabled || recoveryMode)
        problem->restoreReservoirHead(getParam<Scalar>("Restart.ReservoirHeadCm"));

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

    if (recoveryMode) {
        auto fv=localView(*gridGeometry);
        auto vv=localView(gridVariables->curGridVolVars());
        using FS=GetPropType<TypeTag, Properties::FluidSystem>;
        auto check=[](double got,double expected,double atol,double rtol) {
            if(!std::isfinite(got) || std::abs(got-expected)>atol+rtol*std::abs(expected))
                throw std::runtime_error("Recovered primary variables do not reproduce saved VTU fields; stopped before time integration");
        };
        for(const auto& element: elements(leafGridView)) {
            fv.bind(element); vv.bind(element,fv,x);
            for(const auto& scv: scvs(fv)) {
                const auto& v=vv[scv]; const auto& r=recovered[scv.dofIndex()];
                check(v.pressure(0),r.pl,1e-4,1e-8);
                check(v.pressure(1),r.pg,1e-4,1e-8);
                check(v.saturation(1),r.sg,1e-8,1e-6);
                check(v.moleFraction(0,FS::O2Idx),r.o2,1e-13,1e-6);
                check(v.moleFraction(1,FS::O2Idx),r.og,1e-12,1e-6);
                check(v.moleFraction(0,FS::N2Idx),r.n2,1e-13,1e-6);
                check(v.moleFraction(1,FS::N2Idx),r.ng,1e-12,1e-6);
            }
        }
        std::cout << "RECOVERY VALIDATION PASSED: mesh, phase state, pressure, saturation, O2 and N2\n";
    }


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

    // Write the loaded or analytical initial state under the new run name.
    vtkWriter.write(startTime);

    // Read time-loop parameters. DuMuX uses seconds here.
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
            startTime,
            dtInitial,
            tEnd
        );

    timeLoop->setMaxTimeStepSize(maxDt);

    // Let the problem read the currently attempted implicit time level.
    // This remains correct when the Newton solver internally reduces dt and
    // retries the same step.
    problem->setTimeLoop(timeLoop);

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

    // Create the sequential UMFPack direct linear solver.
    using LinearSolver =
        UMFPackIstlSolver<
            LinearSolverTraits<GridGeometry>,
            LinearAlgebraTraitsFromAssembler<Assembler>
        >;

    auto linearSolver =
        std::make_shared<LinearSolver>();

    // Create the nonlinear Newton solver.
    using NonlinearSolver =
        FloodMarChoppedNewtonSolver<Assembler, LinearSolver>;

    NonlinearSolver nonlinearSolver(
        assembler,
        linearSolver
    );

    Scalar nextOutputTime =
        startTime + outputInterval;

    timeLoop->start();

    std::ofstream reservoirHistory(
        problem->name()
        + "_reservoir_history.csv"
    );

    reservoirHistory
        << "time_h,reservoir_head_cm,"
        << "reservoir_storage_m3,"
        << "reservoir_net_flux_m3_s\n";

    reservoirHistory
        << std::setprecision(15);

    std::ofstream injectionHistory(
        problem->name()
        + "_injection_history.csv"
    );

    injectionHistory
        << "time_h,cycle,phase,dt_s,"
        << "applied_boundary_head_cm,"
        << "drywell_net_flux_to_well_m3_s,"
        << "well_to_soil_flux_m3_s,"
        << "soil_to_well_flux_m3_s,"
        << "step_well_to_soil_m3,"
        << "step_soil_to_well_m3,"
        << "step_net_to_soil_m3,"
        << "step_well_storage_change_m3,"
        << "step_external_supply_m3,"
        << "cumulative_wetting_well_to_soil_m3,"
        << "cumulative_wetting_soil_to_well_m3,"
        << "cumulative_wetting_net_to_soil_m3,"
        << "cumulative_drainage_well_to_soil_m3,"
        << "cumulative_drainage_soil_to_well_m3,"
        << "cumulative_drainage_net_to_soil_m3,"
        << "cumulative_total_well_to_soil_m3,"
        << "cumulative_total_soil_to_well_m3,"
        << "cumulative_total_net_to_soil_m3,"
        << "cumulative_external_supply_m3,"
        << "reservoir_state_head_cm,"
        << "reservoir_state_storage_m3\n";

    injectionHistory
        << std::setprecision(15);

    Scalar cumulativeWettingWellToSoilM3 = 0.0;
    Scalar cumulativeWettingSoilToWellM3 = 0.0;
    Scalar cumulativeWettingNetToSoilM3 = 0.0;

    Scalar cumulativeDrainageWellToSoilM3 = 0.0;
    Scalar cumulativeDrainageSoilToWellM3 = 0.0;
    Scalar cumulativeDrainageNetToSoilM3 = 0.0;

    Scalar cumulativeTotalWellToSoilM3 = 0.0;
    Scalar cumulativeTotalSoilToWellM3 = 0.0;
    Scalar cumulativeTotalNetToSoilM3 = 0.0;

    Scalar cumulativeExternalSupplyM3 = 0.0;

    // The production startup ramp begins from zero water head.
    // This variable tracks the physically represented water volume
    // inside the drywell across wetting/drainage transitions.
    Scalar physicalWellStorageM3 = (restartEnabled || recoveryMode) ? problem->reservoirStorageM3() : 0.0;


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

        if (nextOutputTime > timeLoop->time() && nextOutputTime < timeLoop->time()+timeLoop->timeStepSize())
            timeLoop->setTimeStepSize(nextOutputTime-timeLoop->time());
        // Solve the nonlinear two-phase system.
        nonlinearSolver.solve(
            x,
            *timeLoop
        );

        // The Newton solver may have reduced dt internally.
        // Use the actually accepted dt/end time.
        const Scalar acceptedDt =
            timeLoop->timeStepSize();

        const Scalar acceptedEndTime =
            timeLoop->time()
            + acceptedDt;

        // Read-only boundary-flux reconstruction using exactly the
        // hydraulic head that was applied during this accepted step.
        const auto fluxDiag =
            problem->drywellFluxDiagnostics(
                *gridVariables,
                x,
                acceptedEndTime
            );

        const bool externalSupplyActive =
            problem->externalSupplyActive(
                acceptedEndTime
            );

        const int cycle =
            problem->cycleNumber(
                acceptedEndTime
            );

        const Scalar stepWellToSoilM3 =
            fluxDiag.wellToSoilFluxM3s
            *acceptedDt;

        const Scalar stepSoilToWellM3 =
            fluxDiag.soilToWellFluxM3s
            *acceptedDt;

        // Positive means net transfer from the well into soil.
        const Scalar stepNetToSoilM3 =
            -fluxDiag.netFluxToWellM3s
            *acceptedDt;

        Scalar stepWellStorageChangeM3 = 0.0;
        Scalar stepExternalSupplyM3 = 0.0;

        if (externalSupplyActive)
        {
            // During constant-head wetting, external water must provide:
            //   1. net water transferred from well to soil, plus
            //   2. any increase in water stored in the well itself.
            //
            // This captures the initial 0 -> 12 m filling ramp and the
            // refill from residual head back to 12 m at cycles 2 and 3.
            stepWellStorageChangeM3 =
                fluxDiag.boundaryStorageM3
                - physicalWellStorageM3;

            stepExternalSupplyM3 =
                stepNetToSoilM3
                + stepWellStorageChangeM3;

            cumulativeWettingWellToSoilM3 +=
                stepWellToSoilM3;

            cumulativeWettingSoilToWellM3 +=
                stepSoilToWellM3;

            cumulativeWettingNetToSoilM3 +=
                stepNetToSoilM3;

            cumulativeExternalSupplyM3 +=
                stepExternalSupplyM3;
        }
        else
        {
            cumulativeDrainageWellToSoilM3 +=
                stepWellToSoilM3;

            cumulativeDrainageSoilToWellM3 +=
                stepSoilToWellM3;

            cumulativeDrainageNetToSoilM3 +=
                stepNetToSoilM3;
        }

        cumulativeTotalWellToSoilM3 +=
            stepWellToSoilM3;

        cumulativeTotalSoilToWellM3 +=
            stepSoilToWellM3;

        cumulativeTotalNetToSoilM3 +=
            stepNetToSoilM3;

        problem->updateReservoirAfterAcceptedStep(
            *gridVariables,
            x,
            acceptedEndTime,
            acceptedDt
        );

        if (externalSupplyActive)
            physicalWellStorageM3 =
                fluxDiag.boundaryStorageM3;
        else
            physicalWellStorageM3 =
                problem->reservoirStorageM3();

        // Accept the new solution.
        xOld = x;
        gridVariables->advanceTimeStep();

        // Advance simulation time.
        timeLoop->advanceTimeStep();
        timeLoop->reportTimeStep();

        injectionHistory
            << timeLoop->time()/3600.0 << ","
            << cycle << ","
            << (externalSupplyActive
                ? "wetting"
                : "drainage") << ","
            << acceptedDt << ","
            << fluxDiag.boundaryHeadCm << ","
            << fluxDiag.netFluxToWellM3s << ","
            << fluxDiag.wellToSoilFluxM3s << ","
            << fluxDiag.soilToWellFluxM3s << ","
            << stepWellToSoilM3 << ","
            << stepSoilToWellM3 << ","
            << stepNetToSoilM3 << ","
            << stepWellStorageChangeM3 << ","
            << stepExternalSupplyM3 << ","
            << cumulativeWettingWellToSoilM3 << ","
            << cumulativeWettingSoilToWellM3 << ","
            << cumulativeWettingNetToSoilM3 << ","
            << cumulativeDrainageWellToSoilM3 << ","
            << cumulativeDrainageSoilToWellM3 << ","
            << cumulativeDrainageNetToSoilM3 << ","
            << cumulativeTotalWellToSoilM3 << ","
            << cumulativeTotalSoilToWellM3 << ","
            << cumulativeTotalNetToSoilM3 << ","
            << cumulativeExternalSupplyM3 << ","
            << problem->reservoirHeadCm() << ","
            << problem->reservoirStorageM3()
            << "\n";

        if (timeLoop->time() >= 95.0*3600.0)
        {
            reservoirHistory
                << timeLoop->time()/3600.0 << ","
                << problem->reservoirHeadCm() << ","
                << problem->reservoirStorageM3() << ","
                << problem->lastReservoirNetFluxM3s()
                << "\n";

            reservoirHistory.flush();
        }

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


    if (getParam<bool>("Restart.WriteAtEnd",false)) {
        const auto path=getParam<std::string>("Restart.OutputFile");
        std::ofstream snapshot(path);
        snapshot << std::setprecision(17) << "FLOODMAR_RESTART_V1 " << x.size() << " " << timeLoop->time() << "\n";
        for(std::size_t i=0;i<x.size();++i)
            snapshot << i << " " << x[i].state() << " " << x[i][0] << " " << x[i][1] << " " << x[i][2] << "\n";
        snapshot.flush();
        if(!snapshot) throw std::runtime_error("Failed writing restart snapshot");
        std::ofstream headFile(path+".head");
        headFile << std::setprecision(17) << problem->reservoirHeadCm() << "\n";
        headFile.flush();
        if(!headFile) throw std::runtime_error("Failed writing reservoir restart head");
        std::cout << "WROTE RESTART: " << path << " at " << timeLoop->time()/3600.0 << " h\n";
    }

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

