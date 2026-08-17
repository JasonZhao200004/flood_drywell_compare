#ifndef DUMUX_FLOOD_MAR_PROBLEM_HH
#define DUMUX_FLOOD_MAR_PROBLEM_HH

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
    using TimeLoop = Dumux::TimeLoop<Scalar>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

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
    , pondStartupRampTime_(
        getParam<Scalar>("Problem.PondStartupRampTime", 43200.0)
      )
    , pondInjectionHeadCm_(
        getParam<Scalar>("Problem.PondInjectionHead", 5.0)
      )
    , pondDrainageHeadCm_(
        getParam<Scalar>("Problem.PondDrainageHead", -25.0)
      )
    , pondRadius_(
        getParam<Scalar>("Problem.PondRadius", 5.0)
      )
    , pondDepth_(
        getParam<Scalar>("Problem.PondDepth", 0.5)
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
                  << "Pond startup ramp time = "
                  << pondStartupRampTime_ << " s\n"
                  << "Initial O2 profile = "
                  << initialOxygenProfileFile_ << "\n"
                  << "Flood-MAR pond: radius " << pondRadius_
                  << " m, excavation " << pondDepth_ << " m; "
                  << "axisymmetric soil radius 20 m\n"
                  << "HYDRUS-equivalent pond BC for CCTpfa: target-head "
                  << "Darcy flux plus third-type component transport\n"
                  << "Recharge head = " << pondInjectionHeadCm_
                  << " cm; drainage head = " << pondDrainageHeadCm_ << " cm\n"
                  << "Formal pond simulation: 3 cycles, 720 h"
                  << std::endl;
    }

    const std::string& name() const
    { return name_; }

    Scalar temperatureAtPos(const GlobalPosition&) const
    { return 298.15; }

    void setTime(Scalar time)
    { time_ = time; }

    // Keep the boundary synchronized when Newton internally cuts and retries
    // a time step.  Reading t + the currently attempted dt avoids evaluating
    // a retry with the rejected larger time level.
    void setTimeLoop(const std::shared_ptr<TimeLoop>& timeLoop)
    { timeLoop_ = timeLoop; }

    Scalar nextBoundaryEventTime(Scalar currentTimeSeconds) const
    {
        static constexpr std::array<Scalar, 12> eventHours = {
            0.0, 96.0, 101.0, 240.0, 245.0, 336.0,
            341.0, 480.0, 485.0, 576.0, 581.0, 720.0
        };

        Scalar nextEventTime = std::numeric_limits<Scalar>::infinity();

        if (pondStartupRampTime_ > currentTimeSeconds + 1.0e-8)
            nextEventTime = pondStartupRampTime_;

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

        // CCTpfa does not support equation-wise mixed boundary types on one
        // face.  HYDRUS solves the water-flow and solute equations
        // sequentially and can therefore combine VariHead with a third-type
        // solute boundary.  In this fully coupled H2O-N2-O2 system the same
        // physics is written as one pure flux boundary in neumann():
        //   (1) the prescribed pond head defines the exterior pressure;
        //   (2) Darcy's law gives the water flux across the boundary face;
        //   (3) incoming components use the pond-water composition, whereas
        //       outgoing components use the local composition (third type).
        return bcTypes;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (isPondVariHead_(globalPos))
            return pondState_(globalPos);

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
        else if (isPondVariHead_(globalPos))
        {
            // CCTpfa-compatible equivalent of HYDRUS VariHead + third type.
            // The target head is retained exactly as the exterior hydraulic
            // state.  Its pressure difference from the inside cell produces
            // a Darcy water flux.  All component equations are then supplied
            // consistently as flux conditions on this face.
            const int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
            const Scalar density = volVars.density(liquidPhaseIdx);
            const Scalar molarDensity = volVars.molarDensity(liquidPhaseIdx);
            const Scalar mobility = volVars.mobility(liquidPhaseIdx);
            const auto permeability = volVars.permeability();

            constexpr Scalar waterDensity = 1000.0;
            constexpr Scalar gravity = 9.81;
            constexpr Scalar atmosphericPressure = 1.0e5;

            // The pond bottom (r = 0--pondRadius, y = 29.5 m) and submerged
            // part of the vertical rim share the HYDRUS VariHead boundary.
            // The submerged rim height follows the specified positive head,
            // capped by the 0.5 m excavation depth.
            const Scalar localPondHeadCm = effectivePondHeadCm_(globalPos);

            const Scalar boundaryLiquidPressure =
                atmosphericPressure
                + waterDensity*gravity*(localPondHeadCm/100.0);

            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const Scalar distance =
                (insideScv.center() - globalPos).two_norm();
            const auto normal = scvf.unitOuterNormal();

            // Outward liquid Darcy volume flux [m/s]. The gravity vector
            // points in the negative y direction.  effectivePondHeadCm_()
            // already performs the initial ramp, so no second multiplicative
            // ramp is applied to the flux.  The -25 cm stage remains an
            // actual prescribed-head drainage boundary at the pond bottom,
            // matching the HYDRUS VariHead schedule.
            const Scalar liquidVolumeFlux =
                permeability*mobility
                *((volVars.pressure(liquidPhaseIdx)
                   - boundaryLiquidPressure)/distance
                  - density*gravity*normal[1]);

            const Scalar xH2OLocal =
                volVars.moleFraction(liquidPhaseIdx, FluidSystem::H2OIdx);
            const Scalar xN2Local =
                volVars.moleFraction(liquidPhaseIdx, FluidSystem::N2Idx);
            const Scalar xO2Local =
                volVars.moleFraction(liquidPhaseIdx, FluidSystem::O2Idx);

            const Scalar inflowOxygenMolarConcentration =
                pondInfluentOxygenMassConcentration_
                /FluidSystem::molarMass(FluidSystem::O2Idx);

            if (liquidVolumeFlux < 0.0)
            {
                // Inflow: HYDRUS third-type O2 concentration is prescribed.
                // HYDRUS does not solve an N2 solute equation; consequently
                // N2 is kept at the adjacent aqueous composition and H2O is
                // the remainder.  This avoids introducing an artificial N2
                // source solely because DuMuX explicitly represents gas.
                const Scalar xO2In = std::clamp(
                    inflowOxygenMolarConcentration/molarDensity,
                    Scalar(0.0), Scalar(1.0));
                const Scalar xN2In = std::clamp(
                    xN2Local, Scalar(0.0), Scalar(1.0) - xO2In);
                const Scalar xH2OIn = Scalar(1.0) - xN2In - xO2In;

                values[Indices::conti0EqIdx + FluidSystem::H2OIdx] =
                    liquidVolumeFlux*molarDensity*xH2OIn;
                values[Indices::conti0EqIdx + FluidSystem::N2Idx] =
                    liquidVolumeFlux*molarDensity*xN2In;
                values[Indices::conti0EqIdx + FluidSystem::O2Idx] =
                    liquidVolumeFlux*molarDensity*xO2In;
            }
            else
            {
                // Outflow: use the local liquid composition (upwind state).
                values[Indices::conti0EqIdx + FluidSystem::H2OIdx] =
                    liquidVolumeFlux*molarDensity*xH2OLocal;
                values[Indices::conti0EqIdx + FluidSystem::N2Idx] =
                    liquidVolumeFlux*molarDensity*xN2Local;
                values[Indices::conti0EqIdx + FluidSystem::O2Idx] =
                    liquidVolumeFlux*molarDensity*xO2Local;
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
        return stateFromHead_(
            globalPos,
            initialHydraulicHeadCm_(globalPos),
            liquidMoleFractionFromLiquidMassConcentration_(
                initialOxygenMassConcentrationAtElevation_(globalPos[1]))
        );
    }

    void readInitialOxygenProfile_()
    {
        std::ifstream input(initialOxygenProfileFile_);
        if (!input)
            throw std::runtime_error(
                "Cannot open initial oxygen profile: "
                + initialOxygenProfileFile_
            );

        Scalar elevation = 0.0;
        Scalar concentrationGPerCm3 = 0.0;
        while (input >> elevation >> concentrationGPerCm3)
        {
            oxygenProfileElevation_.push_back(elevation);
            // 1 g/cm3 = 1000 kg/m3.
            oxygenProfileMassConcentration_.push_back(
                concentrationGPerCm3*1000.0
            );
        }

        if (oxygenProfileElevation_.size() < 2)
            throw std::runtime_error(
                "Initial oxygen profile must contain at least two rows"
            );

        if (!std::is_sorted(
                oxygenProfileElevation_.begin(),
                oxygenProfileElevation_.end()))
            throw std::runtime_error(
                "Initial oxygen profile elevations must be increasing"
            );
    }

    Scalar initialOxygenMassConcentrationAtElevation_(Scalar elevation) const
    {
        if (elevation <= oxygenProfileElevation_.front())
            return oxygenProfileMassConcentration_.front();
        if (elevation >= oxygenProfileElevation_.back())
            return oxygenProfileMassConcentration_.back();

        const auto upper = std::upper_bound(
            oxygenProfileElevation_.begin(),
            oxygenProfileElevation_.end(),
            elevation
        );
        const std::size_t i1 =
            std::distance(oxygenProfileElevation_.begin(), upper);
        const std::size_t i0 = i1 - 1;
        const Scalar fraction =
            (elevation - oxygenProfileElevation_[i0])
            /(oxygenProfileElevation_[i1]
              - oxygenProfileElevation_[i0]);

        return oxygenProfileMassConcentration_[i0]
            + fraction
              *(oxygenProfileMassConcentration_[i1]
                - oxygenProfileMassConcentration_[i0]);
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

    PrimaryVariables pondState_(const GlobalPosition& globalPos) const
    {
        return stateFromHead_(
            globalPos,
            effectivePondHeadCm_(globalPos),
            liquidMoleFractionFromLiquidMassConcentration_(
                pondInfluentOxygenMassConcentration_
            )
        );
    }

    Scalar initialHydraulicHeadCm_(const GlobalPosition& globalPos) const
    {
        const Scalar depthCm = (30.0 - globalPos[1])*100.0;
        if (depthCm <= 1750.0)
            return -55.0 + depthCm/1750.0*30.0;
        if (depthCm <= 2090.0)
            return -25.0
                   + (depthCm - 1750.0)/(2090.0 - 1750.0)*(-125.0);
        return -150.0
               + (depthCm - 2090.0)/(3000.0 - 2090.0)*50.0;
    }

    Scalar effectivePondHeadCm_(const GlobalPosition& globalPos) const
    {
        const Scalar boundaryTime = boundaryTime_();
        const Scalar target = pondHeadCm_(boundaryTime);
        if (pondStartupRampTime_ <= 0.0
            || boundaryTime >= pondStartupRampTime_)
            return target;

        const Scalar f = std::clamp(
            boundaryTime/pondStartupRampTime_, Scalar(0.0), Scalar(1.0));
        return initialHydraulicHeadCm_(globalPos)
               + f*(target - initialHydraulicHeadCm_(globalPos));
    }

    Scalar pondHeadCm_(Scalar timeSeconds) const
    {
        const Scalar timeHours = timeSeconds/3600.0;
        static constexpr std::array<Scalar, 12> times = {
            0.0, 96.0, 101.0, 240.0, 245.0, 336.0,
            341.0, 480.0, 485.0, 576.0, 581.0, 720.0
        };
        const std::array<Scalar, 12> heads = {
            pondInjectionHeadCm_, pondInjectionHeadCm_,
            pondDrainageHeadCm_, pondDrainageHeadCm_,
            pondInjectionHeadCm_, pondInjectionHeadCm_,
            pondDrainageHeadCm_, pondDrainageHeadCm_,
            pondInjectionHeadCm_, pondInjectionHeadCm_,
            pondDrainageHeadCm_, pondDrainageHeadCm_
        };

        if (timeHours <= times.front()) return heads.front();
        if (timeHours >= times.back()) return heads.back();

        for (std::size_t i = 0; i + 1 < times.size(); ++i)
        {
            if (timeHours >= times[i] && timeHours <= times[i + 1])
            {
                const Scalar f =
                    (timeHours - times[i])/(times[i + 1] - times[i]);
                return heads[i] + f*(heads[i + 1] - heads[i]);
            }
        }
        return heads.back();
    }

    bool isAtmosphere_(const GlobalPosition& globalPos) const
    {
        const bool outsideTop =
            globalPos[1] > 30.0 - eps_
            && globalPos[0] >= pondRadius_ - eps_;
        const bool exposedPondSurface =
            isPondBoundaryCandidate_(globalPos)
            && !isPondVariHead_(globalPos);
        return outsideTop || exposedPondSurface;
    }

    bool isBottom_(const GlobalPosition& globalPos) const
    { return globalPos[1] < eps_; }

    bool isPondVariHead_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const Scalar headCm = pondHeadCm_(boundaryTime_());

        // The pond bottom is VariHead during both recharge (+head) and
        // drainage (-25 cm head).  Only a positive pond depth wets part of
        // the vertical rim.
        const bool pondBottom =
            std::abs(y - pondBottomElevation_) < eps_
            && x <= pondRadius_ + eps_;
        const Scalar wettedRimHeight =
            std::clamp(headCm/100.0, Scalar(0.0), pondDepth_);
        const bool wettedRim =
            std::abs(x - pondRadius_) < eps_
            && y >= pondBottomElevation_ - eps_
            && y <= pondBottomElevation_ + wettedRimHeight + eps_;

        return pondBottom || wettedRim;
    }

    bool isPondBoundaryCandidate_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const bool pondBottom =
            std::abs(y - pondBottomElevation_) < eps_
            && x <= pondRadius_ + eps_;
        const bool completeRim =
            std::abs(x - pondRadius_) < eps_
            && y >= pondBottomElevation_ - eps_
            && y <= pondBottomElevation_ + pondDepth_ + eps_;
        return pondBottom || completeRim;
    }

    static constexpr Scalar eps_ = 1.0e-6;
    static constexpr Scalar pondBottomElevation_ = 29.5;
    Scalar time_ = 0.0;
    std::shared_ptr<const TimeLoop> timeLoop_;
    std::string name_;
    Scalar atmosphericOxygenDiffusionCoefficient_;
    Scalar atmosphericGasPressureFactor_;
    Scalar pondStartupRampTime_;
    Scalar pondInjectionHeadCm_;
    Scalar pondDrainageHeadCm_;
    Scalar pondRadius_;
    Scalar pondDepth_;

    std::string initialOxygenProfileFile_;
    std::vector<Scalar> oxygenProfileElevation_;
    std::vector<Scalar> oxygenProfileMassConcentration_;

    // HYDRUS cValue2: 8e-6 g/cm3 = 0.008 kg/m3 inflowing water.
    static constexpr Scalar pondInfluentOxygenMassConcentration_ = 0.008;
    // 0.000261 g/cm3 = 0.261 kg/m3 gas.
    static constexpr Scalar atmosphericOxygenMassConcentration_ = 0.261;
};

} // namespace Dumux

#endif
