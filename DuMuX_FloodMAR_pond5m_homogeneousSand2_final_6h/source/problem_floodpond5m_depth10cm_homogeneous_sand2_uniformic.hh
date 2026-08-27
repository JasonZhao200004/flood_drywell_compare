#ifndef DUMUX_FLOOD_MAR_PROBLEM_HH_FLOODPOND5M_DEPTH10CM
#define DUMUX_FLOOD_MAR_PROBLEM_HH_FLOODPOND5M_DEPTH10CM

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/binarycoefficients/h2o_o2.hh>

namespace Dumux {

template<class TypeTag>
class WaterAirProblem
: public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum
    {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
        oxygenIdx = Indices::switchIdx + 1
    };

    enum
    {
        wPhaseOnly = Indices::firstPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using TimeLoop = Dumux::TimeLoop<Scalar>;

public:
    WaterAirProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , atmosphericOxygenDiffusionCoefficient_(
        getParam<Scalar>(
            "Problem.AtmosphericOxygenDiffusionCoefficient",
            getParam<Scalar>("Problem.AtmosphericGasTransferCoefficient", 1.0e-6)
        )
      )
    , atmosphericGasPressureFactor_(
        getParam<Scalar>("Problem.AtmosphericGasPressureFactor", 1.0)
      )
    , enableGasPhaseOxygenAdvection_(
        getParam<bool>("Problem.EnableGasPhaseOxygenAdvection", true)
      )
    , drywellStartupRampTime_(
        getParam<Scalar>("Problem.DrywellStartupRampTime", 3600.0)
      )
    , drywellInjectionHeadCm_(
        getParam<Scalar>("Problem.DrywellInjectionHead", 1200.0)
      )
    , initialOxygenProfileFile_(
        getParam<std::string>(
            "Problem.InitialOxygenProfileFile",
            "oxygen_initial_2880.dat"
        )
      )
    {
        FluidSystem::init();
        name_ = getParam<std::string>("Problem.Name", "floodmar");
        readInitialOxygenProfile_();

        std::cout << "Flood-MAR axisymmetric water-N2-O2 model\n"
                  << "Open atmospheric gas-pressure boundary\n"
                  << "O2 surface diffusion coefficient = "
                  << atmosphericOxygenDiffusionCoefficient_ << " m/s\n"
                  << "Gas-pressure transmissibility factor = "
                  << atmosphericGasPressureFactor_ << "\n"
                  << "Drywell startup ramp time = "
                  << drywellStartupRampTime_ << " s\n"
                  << "Drywell injection head = "
                  << drywellInjectionHeadCm_ << " cm\n"
                  << "Drywell boundary mode = legacy whole-face midpoint\n"
                  << "Initial O2 profile = "
                  << initialOxygenProfileFile_ << "\n"
                  << "Formal drywell simulation: 3 cycles, 720 h"
                  << std::endl;
    }

    const std::string& name() const
    { return name_; }

    Scalar temperatureAtPos(const GlobalPosition&) const
    { return 298.15; }

    void setTime(Scalar time)
    { time_ = time; }

    // Keep time-dependent boundary conditions synchronized with the time
    // step currently tried by DuMuX. NewtonSolver::solve(x, timeLoop) may
    // reduce dt internally after a failed Newton attempt. Reading the shared
    // time loop here ensures the retry uses t + the reduced dt, not the end
    // time of the failed larger step.
    void setTimeLoop(const std::shared_ptr<TimeLoop>& timeLoop)
    { timeLoop_ = timeLoop; }

    Scalar nextBoundaryEventTime(Scalar currentTimeSeconds) const
    {
        static constexpr std::array<Scalar, 12> eventHours = {
            0.0, 96.0, 101.0, 240.0, 245.0, 336.0,
            341.0, 480.0, 485.0, 576.0, 581.0, 720.0
        };

        Scalar nextEventTime = std::numeric_limits<Scalar>::infinity();

        // The end of the configurable startup ramp is also a boundary-law
        // change point and should not be crossed by one time step.
        if (drywellStartupRampTime_ > currentTimeSeconds + 1.0e-8)
            nextEventTime = drywellStartupRampTime_;

        for (const Scalar eventHour : eventHours)
        {
            const Scalar eventTime = eventHour*3600.0;
            if (eventTime > currentTimeSeconds + 1.0e-8)
                nextEventTime = std::min(nextEventTime, eventTime);
        }

        return nextEventTime;
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();

        return bcTypes;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (isVariHeadBoundary_(globalPos))
            return drywellState_(globalPos);

        return initial_(globalPos);
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element&,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache&,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];

        if (isAtmosphere_(globalPos))
        {
            const int gasPhaseIdx = FluidSystem::gasPhaseIdx;
            constexpr Scalar atmosphericPressure = 1.0e5;
            constexpr Scalar temperature = 298.15;
            constexpr Scalar gasConstant = 8.31446261815324;
            constexpr Scalar atmosphericWaterVaporPressure = 3169.0;

            const Scalar cO2Atmosphere =
                atmosphericOxygenMassConcentration_
                /FluidSystem::molarMass(FluidSystem::O2Idx);

            const Scalar cH2OAtmosphere =
                atmosphericWaterVaporPressure
                /(gasConstant*temperature);

            const Scalar cDryAtmosphere =
                (atmosphericPressure - atmosphericWaterVaporPressure)
                /(gasConstant*temperature);

            const Scalar cN2Atmosphere =
                std::max(Scalar(0.0), cDryAtmosphere - cO2Atmosphere);

            const Scalar gasMolarDensity = volVars.molarDensity(gasPhaseIdx);
            const Scalar xH2OSoil =
                volVars.moleFraction(gasPhaseIdx, FluidSystem::H2OIdx);
            const Scalar xO2Soil =
                volVars.moleFraction(gasPhaseIdx, FluidSystem::O2Idx);
            const Scalar xN2Soil =
                volVars.moleFraction(gasPhaseIdx, FluidSystem::N2Idx);
            const Scalar cO2Soil = gasMolarDensity*xO2Soil;

            // Outward Darcy gas volume flux [m/s]. Positive values vent
            // soil gas; negative values draw atmospheric gas into the soil.
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const Scalar distance =
                (insideScv.center() - globalPos).two_norm();

            const Scalar gasVolumeFlux =
                atmosphericGasPressureFactor_
                *volVars.permeability()
                *volVars.mobility(gasPhaseIdx)
                *(volVars.pressure(gasPhaseIdx) - atmosphericPressure)
                /distance;

            Scalar advectiveH2OFlux = 0.0;
            Scalar advectiveN2Flux = 0.0;
            Scalar advectiveO2Flux = 0.0;

            if (gasVolumeFlux >= 0.0)
            {
                // Gas leaving the domain carries the local soil composition.
                advectiveH2OFlux =
                    gasVolumeFlux*gasMolarDensity*xH2OSoil;
                advectiveN2Flux =
                    gasVolumeFlux*gasMolarDensity*xN2Soil;
                advectiveO2Flux =
                    gasVolumeFlux*gasMolarDensity*xO2Soil;
            }
            else
            {
                // Gas entering the domain carries the complete atmospheric
                // H2O-N2-O2 composition at 25 degrees Celsius.
                advectiveH2OFlux = gasVolumeFlux*cH2OAtmosphere;
                advectiveN2Flux = gasVolumeFlux*cN2Atmosphere;
                advectiveO2Flux = gasVolumeFlux*cO2Atmosphere;
            }

            // Stagnant-O2 diagnostic: preserve the total advective gas
            // molar flux by treating the suppressed advective O2 as inert N2.
            // Gas pressure can therefore vent normally while O2 itself is
            // transported in the gas phase only by diffusion.
            if (!enableGasPhaseOxygenAdvection_)
            {
                advectiveN2Flux += advectiveO2Flux;
                advectiveO2Flux = 0.0;
            }

            // Molecular O2 exchange with the atmospheric reservoir.
            // Positive is outward according to the DuMuX Neumann convention.
            const Scalar diffusiveO2Flux =
                atmosphericOxygenDiffusionCoefficient_
                *(cO2Soil - cO2Atmosphere);

            values[Indices::conti0EqIdx + FluidSystem::H2OIdx] =
                advectiveH2OFlux;
            values[Indices::conti0EqIdx + FluidSystem::N2Idx] =
                advectiveN2Flux;
            values[Indices::conti0EqIdx + FluidSystem::O2Idx] =
                advectiveO2Flux + diffusiveO2Flux;

            // There is no imposed liquid rainfall or evaporation flux. The
            // H2O flux above is only the water vapor transported with the
            // pressure-driven gas mixture.
        }
        else if (isVariHeadBoundary_(globalPos))
        {
            // CCTpfa requires a pure boundary type on each boundary face.
            // Therefore, impose the variable hydraulic head as a Robin flux
            // for all component balances. This is the flux equivalent of the
            // former pressure/state Dirichlet condition, while retaining the
            // HYDRUS third-type oxygen condition.
            const int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
            const Scalar density = volVars.density(liquidPhaseIdx);
            const Scalar molarDensity = volVars.molarDensity(liquidPhaseIdx);
            const Scalar mobility = volVars.mobility(liquidPhaseIdx);
            const auto permeability = volVars.permeability();

            constexpr Scalar waterDensity = 1000.0;
            constexpr Scalar gravity = 9.81;
            constexpr Scalar atmosphericPressure = 1.0e5;

            const Scalar boundaryTime = boundaryTime_();

            // During startup, ramp the water level/head itself instead of
            // activating the complete drywell wall immediately and merely
            // scaling its flux. This makes the wetted boundary grow upward
            // from the drywell bottom and avoids a simultaneous phase switch
            // along the full 12 m injection interval.
            const Scalar startupFactor =
                drywellStartupRampTime_ > 0.0
                ? std::clamp(boundaryTime/drywellStartupRampTime_, Scalar(0.0), Scalar(1.0))
                : Scalar(1.0);

            const Scalar scheduledDrywellHeadCm = drywellHeadCm_(boundaryTime);
            const Scalar rampedDrywellHeadCm =
                scheduledDrywellHeadCm > 0.0
                ? startupFactor*scheduledDrywellHeadCm
                : scheduledDrywellHeadCm;

            // Legacy boundary used by the earlier HYDRUS-comparison model:
            // evaluate the hydrostatic head at the integration point of the
            // boundary face. A complete face is wet when this midpoint head
            // is positive and dry otherwise. There is deliberately no
            // partial-face wetted-fraction weighting in this variant.
            const Scalar localDrywellHeadCm =
                rampedDrywellHeadCm
                - (globalPos[1] - drywellBottomElevation_)*100.0;

            const Scalar boundaryLiquidPressure =
                atmosphericPressure
                + waterDensity*gravity*(localDrywellHeadCm/100.0);

            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const Scalar distance =
                (insideScv.center() - globalPos).two_norm();
            const auto normal = scvf.unitOuterNormal();

            // HYDRUS option: "Atmospheric BC when the specified nodal
            // pressure head is negative". Rainfall and evaporation are zero,
            // so a face whose midpoint is above the water level receives no
            // imposed liquid-water flux.
            // ------------------------------------------------------------
            // POND_DRYOUT_FLUX_SMOOTHING_0P5CM
            //
            // The CCTpfa implementation represents the prescribed hydraulic
            // head as a Robin/Darcy liquid flux.  For the shallow Flood-MAR
            // pond, the complete horizontal pond bottom is at one elevation,
            // so a hard switch from Robin flux to zero flux at h = 0 causes
            // the whole bottom to deactivate simultaneously.
            //
            // Preserve the original Robin flux for h >= 0.5 cm, but smoothly
            // reduce it to zero over the final 0.5 cm of ponded head.
            //
            // smoothstep(s) = 3*s^2 - 2*s^3
            //
            // h >= 0.5 cm : multiplier = 1
            // 0 < h < 0.5 : multiplier smoothly decreases from 1 to 0
            // h <= 0      : multiplier = 0
            // ------------------------------------------------------------

            const Scalar rawLiquidVolumeFlux =
                permeability*mobility
                    *((volVars.pressure(liquidPhaseIdx)
                       - boundaryLiquidPressure)/distance
                      - density*gravity*normal[1]);

            constexpr Scalar dryOutSmoothHeadCm = 0.5;

            Scalar liquidFluxMultiplier = 0.0;

            if (localDrywellHeadCm >= dryOutSmoothHeadCm)
            {
                liquidFluxMultiplier = 1.0;
            }
            else if (localDrywellHeadCm > 0.0)
            {
                const Scalar s =
                    std::clamp(
                        localDrywellHeadCm/dryOutSmoothHeadCm,
                        Scalar(0.0),
                        Scalar(1.0)
                    );

                liquidFluxMultiplier =
                    s*s*(Scalar(3.0) - Scalar(2.0)*s);
            }

            const Scalar liquidVolumeFlux =
                liquidFluxMultiplier*rawLiquidVolumeFlux;

            std::array<Scalar, FluidSystem::numComponents>
                componentMolarConcentration{};

            if (liquidVolumeFlux >= 0.0)
            {
                // Water leaving the domain carries the complete local liquid
                // composition, including dissolved N2 and O2.
                for (int compIdx = 0;
                     compIdx < FluidSystem::numComponents;
                     ++compIdx)
                {
                    componentMolarConcentration[compIdx] =
                        molarDensity
                        *volVars.moleFraction(liquidPhaseIdx, compIdx);
                }
            }
            else
            {
                // Inflow carries HYDRUS cValue2 = 8e-6 g/cm3
                // = 0.008 kg/m3. Incoming water contains no imposed N2;
                // the remaining liquid molar concentration is H2O.
                const Scalar oxygenMolarConcentration =
                    drywellInfluentOxygenMassConcentration_
                    /FluidSystem::molarMass(FluidSystem::O2Idx);

                componentMolarConcentration[FluidSystem::O2Idx] =
                    oxygenMolarConcentration;
                componentMolarConcentration[FluidSystem::N2Idx] = 0.0;
                componentMolarConcentration[FluidSystem::H2OIdx] =
                    std::max(
                        Scalar(0.0),
                        molarDensity - oxygenMolarConcentration
                    );
            }

            for (int compIdx = 0;
                 compIdx < FluidSystem::numComponents;
                 ++compIdx)
            {
                values[Indices::conti0EqIdx + compIdx] =
                    liquidVolumeFlux
                    *componentMolarConcentration[compIdx];
            }
        }
        else if (isBottom_(globalPos))
        {
            constexpr Scalar gravity = 9.81;
            const int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;

            const Scalar density = volVars.density(liquidPhaseIdx);
            const Scalar molarDensity = volVars.molarDensity(liquidPhaseIdx);
            const Scalar mobility = volVars.mobility(liquidPhaseIdx);
            const auto permeability = volVars.permeability();

            // Free drainage: unit hydraulic gradient.
            const Scalar liquidVolumeFlux =
                permeability*mobility*density*gravity;

            for (int compIdx = 0;
                 compIdx < FluidSystem::numComponents;
                 ++compIdx)
            {
                values[Indices::conti0EqIdx + compIdx] =
                    liquidVolumeFlux
                    *molarDensity
                    *volVars.moleFraction(liquidPhaseIdx, compIdx);
            }
        }

        return values;
    }

    // Allow the prescribed drywell liquid flux to participate in
    // near-boundary mechanical-dispersion velocity reconstruction.
    // This does not add a separate boundary dispersive component flux.
    bool isDrywellLiquidBoundaryForDispersion(
        const GlobalPosition& globalPos) const
    { return isVariHeadBoundary_(globalPos); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return initial_(globalPos); }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element&,
                       const FVElementGeometry&,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector values(0.0);
        const auto& volVars = elemVolVars[scv];
        const int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
        const int gasPhaseIdx = FluidSystem::gasPhaseIdx;

        const Scalar porosity = volVars.porosity();
        const Scalar liquidSaturation = volVars.saturation(liquidPhaseIdx);
        const Scalar gasSaturation = volVars.saturation(gasPhaseIdx);
        const Scalar waterContent = porosity*liquidSaturation;
        const Scalar moistureFactor = moistureDecayFactor_(waterContent);

        constexpr Scalar oxygenRegularization = 1.0e-12;
        const auto regularizedOxygenMoleFraction =
            [=](Scalar x)
            {
                const Scalar ratio = x/oxygenRegularization;
                return x*(Scalar(1.0) - std::exp(-ratio*ratio));
            };

        const Scalar xO2Liquid = regularizedOxygenMoleFraction(
            volVars.moleFraction(liquidPhaseIdx, FluidSystem::O2Idx));
        const Scalar xO2Gas = regularizedOxygenMoleFraction(
            volVars.moleFraction(gasPhaseIdx, FluidSystem::O2Idx));

        const Scalar liquidOxygenInventory =
            porosity*liquidSaturation
            *volVars.molarDensity(liquidPhaseIdx)*xO2Liquid;
        const Scalar gasOxygenInventory =
            porosity*gasSaturation
            *volVars.molarDensity(gasPhaseIdx)*xO2Gas;

        // HYDRUS first-order constants, converted from 1/h to 1/s.
        constexpr Scalar liquidDecayRate = 0.06/3600.0;
        constexpr Scalar gasDecayRate = 0.002/3600.0;

        values[Indices::conti0EqIdx + FluidSystem::O2Idx] =
            -moistureFactor
            *(liquidDecayRate*liquidOxygenInventory
              + gasDecayRate*gasOxygenInventory);

        return values;
    }

private:
    Scalar boundaryTime_() const
    {
        if (timeLoop_)
            return timeLoop_->time() + timeLoop_->timeStepSize();

        return time_;
    }

    static Scalar liquidMoleFractionFromLiquidMassConcentration_(Scalar c)
    {
        constexpr Scalar waterDensity = 1000.0;
        const Scalar oxygenMolarConcentration =
            c/FluidSystem::molarMass(FluidSystem::O2Idx);
        const Scalar waterMolarConcentration =
            waterDensity/FluidSystem::molarMass(FluidSystem::H2OIdx);
        return oxygenMolarConcentration
               /(waterMolarConcentration + oxygenMolarConcentration);
    }

    static Scalar moistureDecayFactor_(Scalar theta)
    {
        constexpr Scalar minimumDryRate = 0.1;
        constexpr Scalar theta0 = 0.06;
        constexpr Scalar theta1 = 0.15;
        constexpr Scalar theta2 = 0.25;
        constexpr Scalar theta3 = 0.37;
        constexpr Scalar minimumWetRate = 0.3;

        if (theta <= theta0) return minimumDryRate;
        if (theta < theta1)
        {
            const Scalar f = (theta - theta0)/(theta1 - theta0);
            return minimumDryRate + f*(1.0 - minimumDryRate);
        }
        if (theta <= theta2) return 1.0;
        if (theta < theta3)
        {
            const Scalar f = (theta - theta2)/(theta3 - theta2);
            return 1.0 + f*(minimumWetRate - 1.0);
        }
        return minimumWetRate;
    }

    PrimaryVariables initial_(const GlobalPosition& globalPos) const
    {
        // HOMOGENEOUS SAND2 WATER IC
        // theta = 0.11692665151505993; equivalent Sand2 head = -29.770236344729383 cm.
        const Scalar headCm = -29.770236344729383;

        return stateFromHead_(
            globalPos,
            headCm,
            initialOxygenMoleFractionAtPosition_(globalPos)
        );
    }

    void readInitialOxygenProfile_()
    {
        std::ifstream input(initialOxygenProfileFile_);
        if (!input)
            throw std::runtime_error(
                "Cannot open 2-D initial oxygen field: "
                + initialOxygenProfileFile_
            );

        Scalar radius = 0.0;
        Scalar elevation = 0.0;
        Scalar liquidOxygenMoleFraction = 0.0;
        std::string line;

        while (std::getline(input, line))
        {
            if (line.empty() || line[0] == '#')
                continue;

            std::istringstream stream(line);
            if (!(stream >> radius >> elevation >> liquidOxygenMoleFraction))
                continue;

            oxygenProfileRadius_.push_back(radius);
            oxygenProfileElevation_.push_back(elevation);
            oxygenProfileLiquidMoleFraction_.push_back(liquidOxygenMoleFraction);
        }

        if (oxygenProfileRadius_.empty())
            throw std::runtime_error("2-D initial oxygen field contains no usable rows");

        std::cout << "Loaded 2-D O2 IC cells = "
                  << oxygenProfileRadius_.size() << std::endl;
    }

    Scalar initialOxygenMoleFractionAtPosition_(const GlobalPosition& globalPos) const
    {
        const Scalar radius = globalPos[0];
        const Scalar elevation = globalPos[1];
        Scalar bestDistanceSquared = std::numeric_limits<Scalar>::max();
        Scalar bestValue = oxygenProfileLiquidMoleFraction_.front();

        for (std::size_t i = 0; i < oxygenProfileRadius_.size(); ++i)
        {
            const Scalar dr = radius - oxygenProfileRadius_[i];
            const Scalar dz = elevation - oxygenProfileElevation_[i];
            const Scalar d2 = dr*dr + dz*dz;
            if (d2 < bestDistanceSquared)
            {
                bestDistanceSquared = d2;
                bestValue = oxygenProfileLiquidMoleFraction_[i];
                if (d2 < 1.0e-24) break;
            }
        }
        return bestValue;
    }

    PrimaryVariables stateFromHead_(const GlobalPosition& globalPos,
                                    Scalar headCm,
                                    Scalar liquidOxygenMoleFraction) const
    {
        PrimaryVariables priVars(0.0);
        constexpr Scalar atmosphericPressure = 1.0e5;
        constexpr Scalar waterDensity = 1000.0;
        constexpr Scalar gravity = 9.81;

        const Scalar liquidPressure =
            atmosphericPressure
            + waterDensity*gravity*(headCm/100.0);

        if (headCm >= 0.0)
        {
            priVars.setState(wPhaseOnly);
            priVars[pressureIdx] = liquidPressure;
            priVars[switchIdx] = 0.0;
            priVars[oxygenIdx] = liquidOxygenMoleFraction;
        }
        else
        {
            const Scalar capillaryPressure =
                atmosphericPressure - liquidPressure;
            const auto interaction =
                this->spatialParams().fluidMatrixInteractionAtPos(globalPos);
            Scalar liquidSaturation = interaction.sw(capillaryPressure);
            liquidSaturation = std::clamp(
                liquidSaturation, Scalar(0.0), Scalar(1.0));

            priVars.setState(bothPhases);
            priVars[pressureIdx] = liquidPressure;
            priVars[switchIdx] = 1.0 - liquidSaturation;
            priVars[oxygenIdx] = liquidOxygenMoleFraction;
        }

        return priVars;
    }

    PrimaryVariables drywellState_(const GlobalPosition& globalPos) const
    {
        return stateFromHead_(
            globalPos,
            drywellHeadCm_(boundaryTime_()),
            initialOxygenMoleFractionAtPosition_(globalPos)
        );
    }

    Scalar drywellHeadCm_(Scalar timeSeconds) const
    {
        const Scalar timeHours =
            timeSeconds/3600.0;

        // ------------------------------------------------------------
        // FLOOD-MAR shallow-pond schedule
        //
        // The drywell baseline used a 5 h transition directly between
        // positive ponding head and -25 cm:
        //
        //     +5 cm -> -25 cm
        //
        // For the 5-cm-deep Flood-MAR pond this caused the complete
        // horizontal pond bottom to cross h = 0 very early in the
        // transition and triggered severe nonlinear convergence
        // problems.
        //
        // Here the physical ponded-water stage is allowed to drain
        // smoothly from +5 cm to 0 cm over the complete 5 h transition.
        //
        // After the pond has reached zero water depth, the scheduled
        // dry-period head is again -25 cm.
        //
        // Filling transitions likewise begin hydraulically at 0 cm and
        // rise smoothly to +5 cm.
        //
        // All other model physics and numerical settings are unchanged.
        // ------------------------------------------------------------

        const Scalar injectionHeadCm =
            drywellInjectionHeadCm_;

        constexpr Scalar dryHeadCm = -25.0;

        const auto smoothStep =
            [](Scalar u)
            {
                u = std::clamp(
                    u,
                    Scalar(0.0),
                    Scalar(1.0)
                );

                return
                    u*u
                    *(Scalar(3.0) - Scalar(2.0)*u);
            };

        const auto smoothDown =
            [&](Scalar startHour)
            {
                const Scalar u =
                    (timeHours - startHour)/5.0;

                return
                    injectionHeadCm
                    *(Scalar(1.0) - smoothStep(u));
            };

        const auto smoothUp =
            [&](Scalar startHour)
            {
                const Scalar u =
                    (timeHours - startHour)/5.0;

                return
                    injectionHeadCm
                    *smoothStep(u);
            };

        // Cycle 1
        if (timeHours <= 96.0)
            return injectionHeadCm;

        if (timeHours <= 101.0)
            return smoothDown(96.0);

        if (timeHours < 240.0)
            return dryHeadCm;

        // Cycle 2 filling
        if (timeHours <= 245.0)
            return smoothUp(240.0);

        if (timeHours <= 336.0)
            return injectionHeadCm;

        if (timeHours <= 341.0)
            return smoothDown(336.0);

        if (timeHours < 480.0)
            return dryHeadCm;

        // Cycle 3 filling
        if (timeHours <= 485.0)
            return smoothUp(480.0);

        if (timeHours <= 576.0)
            return injectionHeadCm;

        if (timeHours <= 581.0)
            return smoothDown(576.0);

        return dryHeadCm;
    }

    bool isAtmosphere_(const GlobalPosition& globalPos) const
    { return globalPos[1] > 30.0 - eps_; }

    bool isBottom_(const GlobalPosition& globalPos) const
    { return globalPos[1] < eps_; }

    bool isVariHeadBoundary_(
        const GlobalPosition& globalPos
    ) const
    {
        const Scalar r = globalPos[0];
        const Scalar z = globalPos[1];

        // ------------------------------------------------------------
        // Flood-MAR pond geometry
        //
        // radius       = 5.0 m
        // excavation   = 0.10 m
        // pond bottom  = z 29.90 m
        // head         = 5 cm
        // water level  = z 29.95 m
        //
        // VariHead:
        //   1. complete horizontal pond bottom, r = 0..5 m
        //   2. ONLY lower 5 cm of vertical pond wall,
        //      r = 5 m, z = 29.90..29.95 m
        //
        // Upper 5 cm of pond wall is intentionally NOT VariHead
        // and therefore remains zero-Neumann/no-flow, matching
        // the requested boundary definition.
        // ------------------------------------------------------------

        const bool pondBottom =
            std::abs(z - 29.90) < eps_
            && r <= 5.0 + eps_;

        const bool wettedPondSide =
            std::abs(r - 5.0) < eps_
            && z >= 29.90 - eps_
            && z <= 29.95 + eps_;

        return pondBottom || wettedPondSide;
    }


    static constexpr Scalar eps_ = 1.0e-6;
    static constexpr Scalar drywellRadius_ = 5.0;
    static constexpr Scalar drywellBottomElevation_ = 29.90;
    Scalar time_ = 0.0;
    std::shared_ptr<const TimeLoop> timeLoop_;
    std::string name_;
    Scalar atmosphericOxygenDiffusionCoefficient_;
    Scalar atmosphericGasPressureFactor_;
    bool enableGasPhaseOxygenAdvection_;
    Scalar drywellStartupRampTime_;
    Scalar drywellInjectionHeadCm_;

    std::string initialOxygenProfileFile_;
    std::vector<Scalar> oxygenProfileRadius_;
    std::vector<Scalar> oxygenProfileElevation_;
    std::vector<Scalar> oxygenProfileLiquidMoleFraction_;

    // HYDRUS cValue2: 8e-6 g/cm3 = 0.008 kg/m3 inflowing water.
    static constexpr Scalar drywellInfluentOxygenMassConcentration_ = 0.008;
    // 0.000261 g/cm3 = 0.261 kg/m3 gas.
    static constexpr Scalar atmosphericOxygenMassConcentration_ = 0.261;
};

} // namespace Dumux

#endif
