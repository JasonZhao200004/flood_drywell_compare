#ifndef DUMUX_FLOOD_MAR_OXYGEN_INITIAL_PROPERTIES_HH
#define DUMUX_FLOOD_MAR_OXYGEN_INITIAL_PROPERTIES_HH

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>

#include "dispersionmodeltraits.hh"
#include "transientcctpfadispersionflux.hh"
#include "floodmarprimaryvariableswitch.hh"
#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {

struct FloodMar
{
    using InheritsFrom =
        std::tuple<TwoPNC, CCTpfaModel>;
};

} // namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::FloodMar>
{
    using type = Dune::UGGrid<2>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::FloodMar>
{
    using type = WaterAirProblem<TypeTag>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FloodMar>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = FluidSystems::H2ON2O2<Scalar>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FloodMar>
{
private:
    using GridGeometry =
        GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = WaterAirSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::FloodMar>
{
    static constexpr bool value = true;
};

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::FloodMar>
{
private:
    using BaseTraits =
        GetPropType<TypeTag, Properties::BaseModelTraits>;

public:
    using type = FloodMarDispersionModelTraits<BaseTraits>;
};

template<class TypeTag>
struct DispersionFluxType<TypeTag, TTag::FloodMar>
{
    using type = FloodMarTransientCCTpfaDispersionFlux<
        TypeTag,
        ReferenceSystemFormulation::molarAveraged
    >;
};

template<class TypeTag>
struct Formulation<TypeTag, TTag::FloodMar>
{
    static constexpr auto value = TwoPFormulation::p0s1;
};

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FloodMar>
{
private:
    static constexpr bool enableCache =
        getPropValue<TypeTag, Properties::EnableGridGeometryCache>();

    using GridView =
        typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    struct GGTraits
    : public CCTpfaDefaultGridGeometryTraits<GridView>
    {
        using Extrusion = RotationalExtrusion<0>;
    };

public:
    using type = CCTpfaFVGridGeometry<
        GridView,
        enableCache,
        GGTraits
    >;
};

/*
 * Use the normal 2pnc volume-variable implementation, but export the
 * project-local primary-variable switch.  Reconstructing the traits here
 * mirrors the standard TwoPNC VolumeVariables property.
 */
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::FloodMar>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<
        TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DM = typename GetPropType<
        TypeTag, Properties::GridGeometry>::DiscretizationMethod;

    static constexpr bool enableIS = getPropValue<
        TypeTag, Properties::EnableBoxInterfaceSolver>();

    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<
        PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<
        TypeTag, Properties::EffectiveDiffusivityModel>;

    template<class BaseTraitsT, class DTT, class EDMT>
    struct NCTraits : public BaseTraitsT
    {
        using DiffusionType = DTT;
        using EffectiveDiffusivityModel = EDMT;
    };

public:
    using type = FloodMarTwoPNCVolumeVariables<
        NCTraits<BaseTraits, DT, EDM>>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FloodMar>
{
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FloodMar>
{
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FloodMar>
{
    static constexpr bool value = true;
};

} // namespace Dumux::Properties

#endif
