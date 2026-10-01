#ifndef DUMUX_FLOOD_MAR_SPATIAL_PARAMS_HH
#define DUMUX_FLOOD_MAR_SPATIAL_PARAMS_HH

#include <array>
#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class WaterAirSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<
      GridGeometry,
      Scalar,
      WaterAirSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = WaterAirSpatialParams<GridGeometry, Scalar>;
    using ParentType =
        FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    // Scheidegger mechanical dispersion parameters [m].
    // These correspond to the HYDRUS longitudinal/transverse
    // dispersivities used by the Flood-MAR model.
    std::array<Scalar, 2> dispersionAlphas(
        const GlobalPosition&,
        const int,
        const int
    ) const
    {
        return {
            getParam<Scalar>(
                "Dispersion.LongitudinalDispersivity",
                0.50
            ),
            getParam<Scalar>(
                "Dispersion.TransverseDispersivity",
                0.05
            )
        };
    }

    // Overload required by the thermal-dispersion interface.
    std::array<Scalar, 2> dispersionAlphas(
        const GlobalPosition&,
        const int
    ) const
    {
        return {
            getParam<Scalar>(
                "Dispersion.LongitudinalDispersivity",
                0.50
            ),
            getParam<Scalar>(
                "Dispersion.TransverseDispersivity",
                0.05
            )
        };
    }


    WaterAirSpatialParams(
        std::shared_ptr<const GridGeometry> gridGeometry
    )
    : ParentType(gridGeometry)
    , pcKrSwCurves_{
        PcKrSwCurve("SpatialParams.Sand1"),
        PcKrSwCurve("SpatialParams.LoamySand"),
        PcKrSwCurve("SpatialParams.Sand2"),
        PcKrSwCurve("SpatialParams.Loam"),
        PcKrSwCurve("SpatialParams.SandyLoam"),
        PcKrSwCurve("SpatialParams.Sand3")
      }
    {
        // Intrinsic permeability [m2], converted from HYDRUS Ks.
        permeability_ = {
            1.127e-11, // Sand 1:      Ks = 39.8 cm/h
            3.709e-12, // Loamy Sand:  Ks = 13.1 cm/h
            1.925e-12, // Sand 2:      Ks = 6.8 cm/h
            1.161e-13, // Loam:        Ks = 0.41 cm/h
            2.633e-13, // Sandy Loam:  Ks = 0.93 cm/h
            1.812e-12  // Sand 3:      Ks = 6.4 cm/h
        };

        // HYDRUS saturated water contents are used as porosities.
        porosity_ = {
            0.380733,
            0.384324,
            0.381022,
            0.395717,
            0.389719,
            0.384324
        };
    }

    template<class ElementSolution>
    PermeabilityType permeability(
        const Element& element,
        const SubControlVolume& scv,
        const ElementSolution& elemSol
    ) const
    {
        return permeability_[materialIndex_(
            element.geometry().center()
        )];
    }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_[materialIndex_(globalPos)];
    }

    template<class ElementSolution>
    auto fluidMatrixInteraction(
        const Element& element,
        const SubControlVolume& scv,
        const ElementSolution& elemSol
    ) const
    {
        const auto material = materialIndex_(
            element.geometry().center()
        );

        return makeFluidMatrixInteraction(
            pcKrSwCurves_[material]
        );
    }

    auto fluidMatrixInteractionAtPos(
        const GlobalPosition& globalPos
    ) const
    {
        return makeFluidMatrixInteraction(
            pcKrSwCurves_[materialIndex_(globalPos)]
        );
    }

    template<class FluidSystem>
    int wettingPhaseAtPos(
        const GlobalPosition& globalPos
    ) const
    {
        return FluidSystem::phase0Idx;
    }

    // Kept because the original problem calls this function.
    void plotMaterialLaw() {}

private:
    std::size_t materialIndex_(
        const GlobalPosition& globalPos
    ) const
    {
        const Scalar y = globalPos[1];

        // y is elevation above the 30 m model bottom.
        if (y >= 22.0)
            return 0; // Sand 1: 0 to -8.0 m

        if (y >= 18.5)
            return 1; // Loamy Sand: -8.0 to -11.5 m

        if (y >= 12.5)
            return 2; // Sand 2: -11.5 to -17.5 m

        if (y >= 8.5)
            return 3; // Loam: -17.5 to -21.5 m

        if (y >= 4.0)
            return 4; // Sandy Loam: -21.5 to -26.0 m

        return 5; // Sand 3: -26.0 to -30.0 m
    }

    std::array<Scalar, 6> permeability_;
    std::array<Scalar, 6> porosity_;
    std::array<PcKrSwCurve, 6> pcKrSwCurves_;
};

} // namespace Dumux

#endif