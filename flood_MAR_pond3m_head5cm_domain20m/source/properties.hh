#ifndef DUMUX_FLOOD_MAR_POND_PROPERTIES_HH
#define DUMUX_FLOOD_MAR_POND_PROPERTIES_HH

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>

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

// Two-dimensional unstructured triangular grid
template<class TypeTag>
struct Grid<TypeTag, TTag::FloodMar>
{
    using type = Dune::UGGrid<2>;
};

// Flood-MAR pond recharge problem
template<class TypeTag>
struct Problem<TypeTag, TTag::FloodMar>
{
    using type = WaterAirProblem<TypeTag>;
};

// Water-nitrogen-oxygen fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FloodMar>
{
private:
    using Scalar =
        GetPropType<TypeTag, Properties::Scalar>;

public:
    using type =
        FluidSystems::H2ON2O2<Scalar>;
};

// Six-layer soil spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FloodMar>
{
private:
    using GridGeometry =
        GetPropType<TypeTag, Properties::GridGeometry>;

    using Scalar =
        GetPropType<TypeTag, Properties::Scalar>;

public:
    using type =
        WaterAirSpatialParams<GridGeometry, Scalar>;
};

// Use mole fractions for component primary variables
template<class TypeTag>
struct UseMoles<TypeTag, TTag::FloodMar>
{
    static constexpr bool value = true;
};

// Liquid pressure and gas saturation formulation:
// primary pressure   = liquid pressure
// saturation variable = gas saturation
template<class TypeTag>
struct Formulation<TypeTag, TTag::FloodMar>
{
    static constexpr auto value =
        TwoPFormulation::p0s1;
};

// Axisymmetric geometry:
// x is radial distance;
// rotation is around the y-axis.
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FloodMar>
{
private:
    static constexpr bool enableCache =
        getPropValue<
            TypeTag,
            Properties::EnableGridGeometryCache
        >();

    using GridView =
        typename GetPropType<
            TypeTag,
            Properties::Grid
        >::LeafGridView;

    struct GGTraits
    : public CCTpfaDefaultGridGeometryTraits<GridView>
    {
        using Extrusion =
            RotationalExtrusion<0>;
    };

public:
    using type =
        CCTpfaFVGridGeometry<
            GridView,
            enableCache,
            GGTraits
        >;
};

// Enable geometry caching
template<class TypeTag>
struct EnableGridGeometryCache<
    TypeTag,
    TTag::FloodMar
>
{
    static constexpr bool value = true;
};

// Enable volume-variable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<
    TypeTag,
    TTag::FloodMar
>
{
    static constexpr bool value = true;
};

// Enable flux-variable caching
template<class TypeTag>
struct EnableGridFluxVariablesCache<
    TypeTag,
    TTag::FloodMar
>
{
    static constexpr bool value = true;
};

} // namespace Dumux::Properties

#endif
