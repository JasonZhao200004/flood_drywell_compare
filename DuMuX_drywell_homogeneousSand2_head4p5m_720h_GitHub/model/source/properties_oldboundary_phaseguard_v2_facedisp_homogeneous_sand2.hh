#ifndef DUMUX_FLOODMAR_PHASEGUARD_V2_FACE_DISP_PROPERTIES_HH_HOMOGENEOUS_SAND2
#define DUMUX_FLOODMAR_PHASEGUARD_V2_FACE_DISP_PROPERTIES_HH_HOMOGENEOUS_SAND2

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>

#include "dispersionmodeltraits.hh"
#include "transientcctpfadispersionflux_facesymmetric_v2.hh"
#include "floodmarprimaryvariableswitch_phaseguard_v2.hh"
#include "spatialparams_homogeneous_sand2.hh"
#include "problem_oldboundary_homogeneous_sand2_uniformic.hh"

namespace Dumux::Properties {

namespace TTag {

struct FloodMarOldBoundaryPhaseGuardV2FaceDisp
{
    using InheritsFrom =
        std::tuple<TwoPNC, CCTpfaModel>;
};

} // namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
    using type = Dune::UGGrid<2>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
    using type = WaterAirProblem<TypeTag>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = FluidSystems::H2ON2O2<Scalar>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
private:
    using GridGeometry =
        GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = WaterAirSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
    static constexpr bool value = true;
};

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
private:
    using BaseTraits =
        GetPropType<TypeTag, Properties::BaseModelTraits>;

public:
    using type = FloodMarDispersionModelTraits<BaseTraits>;
};

template<class TypeTag>
struct DispersionFluxType<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
    using type = FloodMarTransientCCTpfaDispersionFlux<
        TypeTag,
        ReferenceSystemFormulation::molarAveraged
    >;
};

template<class TypeTag>
struct Formulation<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
    static constexpr auto value = TwoPFormulation::p0s1;
};

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
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
struct VolumeVariables<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
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
    using type = FloodMarPhaseGuardV2TwoPNCVolumeVariables<
        NCTraits<BaseTraits, DT, EDM>>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FloodMarOldBoundaryPhaseGuardV2FaceDisp>
{
    static constexpr bool value = true;
};

} // namespace Dumux::Properties

#endif
