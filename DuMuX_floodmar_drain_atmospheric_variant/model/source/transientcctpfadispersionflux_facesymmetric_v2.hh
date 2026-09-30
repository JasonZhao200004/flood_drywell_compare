// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef DUMUX_FLOODMAR_TRANSIENT_CCTPFA_FACE_SYMMETRIC_V2_HH
#define DUMUX_FLOODMAR_TRANSIENT_CCTPFA_FACE_SYMMETRIC_V2_HH

#include <algorithm>
#include <cmath>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

/*
 * Project-local transient, two-phase mechanical-dispersion flux for CCTpfa.
 *
 * The stock DuMuX version shipped with this project permits CCTpfa dispersion
 * only for stationary velocity fields and its Scheidegger implementation only
 * for one phase.  This class deliberately avoids changing the DuMuX sources.
 *
 * A cell velocity is reconstructed by a least-squares fit of the phase pore-
 * velocity normal components on all faces of the current cell.  The
 * Scheidegger tensor is then
 *
 *   D = alpha_T |v| I + (alpha_L-alpha_T) v v^T/|v|.
 *
 * By default mechanical dispersion is applied only to the liquid phase, which
 * matches the HYDRUS solute-transport dispersivity parameters.  Gas transport
 * retains its molecular diffusion.  Set Dispersion.ApplyToGasPhase=true only
 * for a deliberate sensitivity test.
 */
template<class TypeTag,
         ReferenceSystemFormulation referenceSystem =
             ReferenceSystemFormulation::molarAveraged>
class FloodMarTransientCCTpfaDispersionFlux
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ElementVolumeVariables =
        typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridFluxVariablesCache =
        GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using ElementFluxVariablesCache =
        typename GridFluxVariablesCache::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numComponents = ModelTraits::numFluidComponents();

    using Matrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using Vector = Dune::FieldVector<Scalar, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

public:
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    static ComponentFluxVector compositionalDispersionFlux(
        const Problem& problem,
        const Element& element,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf,
        const int phaseIdx,
        const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        ComponentFluxVector componentFlux(0.0);

        // All external component fluxes (including the drywell Robin flux)
        // are prescribed by problem.hh.  Do not add an independent
        // mechanical-dispersion flux on those same boundary faces.
        if (scvf.boundary())
            return componentFlux;

        static const bool applyToGas =
            getParam<bool>("Dispersion.ApplyToGasPhase", false);
        if (phaseIdx == FluidSystem::gasPhaseIdx && !applyToGas)
            return componentFlux;

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const Scalar rho = 0.5*(density_(insideVolVars, phaseIdx)
                              + density_(outsideVolVars, phaseIdx));

        const Matrix dispersionTensor = dispersionTensor_(
            problem, element, fvGeometry, elemVolVars,
            scvf, elemFluxVarsCache, phaseIdx
        );
        const Scalar tij = transmissibility_(
            fvGeometry, elemVolVars, scvf, dispersionTensor
        );

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

            const Scalar xInside = fraction_(insideVolVars, phaseIdx, compIdx);
            const Scalar xOutside = fraction_(outsideVolVars, phaseIdx, compIdx);
            componentFlux[compIdx] = rho*tij*(xInside - xOutside);

            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx))
                    componentFlux[FluidSystem::getMainComponent(phaseIdx)]
                        -= componentFlux[compIdx];
        }

        return componentFlux;
    }

    static Scalar thermalDispersionFlux(
        const Problem&, const Element&, const FVElementGeometry&,
        const ElementVolumeVariables&, const SubControlVolumeFace&,
        const int, const ElementFluxVariablesCache&)
    { return 0.0; }

private:
    template<class VolumeVariables>
    static Scalar density_(const VolumeVariables& volVars, const int phaseIdx)
    {
        if constexpr (referenceSystem == ReferenceSystemFormulation::molarAveraged)
            return volVars.molarDensity(phaseIdx);
        else
            return volVars.density(phaseIdx);
    }

    template<class VolumeVariables>
    static Scalar fraction_(const VolumeVariables& volVars,
                            const int phaseIdx,
                            const int compIdx)
    {
        if constexpr (referenceSystem == ReferenceSystemFormulation::molarAveraged)
            return volVars.moleFraction(phaseIdx, compIdx);
        else
            return volVars.massFraction(phaseIdx, compIdx);
    }

    static Matrix dispersionTensor_(
        const Problem& problem,
        const Element& element,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf,
        const ElementFluxVariablesCache& elemFluxVarsCache,
        const int phaseIdx)
    {
        static_assert(
            dimWorld == 2,
            "Face-symmetric FloodMAR dispersion implementation is 2-D."
        );

        Matrix tensor(0.0);

        // Mechanical dispersion remains internal-face only.
        if (scvf.boundary())
            return tensor;

        const Scalar area =
            Extrusion::area(fvGeometry, scvf);

        if (!std::isfinite(area)
            || area <= Scalar(1.0e-30))
            return tensor;

        const auto insideIdx = scvf.insideScvIdx();
        const auto outsideIdx = scvf.outsideScvIdx();

        const auto& insideScv =
            fvGeometry.scv(insideIdx);

        const auto& outsideScv =
            fvGeometry.scv(outsideIdx);

        const auto& insideVolVars =
            elemVolVars[insideIdx];

        const auto& outsideVolVars =
            elemVolVars[outsideIdx];

        // ----------------------------------------------------------
        // Exact CCTpfa Darcy normal flux q_n [m/s].
        // ----------------------------------------------------------
        const Scalar mobility =
            Scalar(0.5)
            *(insideVolVars.mobility(phaseIdx)
              + outsideVolVars.mobility(phaseIdx));

        const Scalar qNormal =
            AdvectionType::flux(
                problem,
                element,
                fvGeometry,
                elemVolVars,
                scvf,
                phaseIdx,
                elemFluxVarsCache
            )
            *mobility/area;

        if (!std::isfinite(qNormal))
            return tensor;

        // ----------------------------------------------------------
        // Construct ONE face-symmetric full Darcy velocity.
        //
        // Direction: line connecting neighboring cell centers.
        //
        // Scale it such that
        //
        //        velocity dot n = qNormal
        //
        // This keeps the exact hydraulic face flux while avoiding the
        // previous element-dependent least-squares tensor.
        // ----------------------------------------------------------
        Vector direction(0.0);

        for (int i = 0; i < dimWorld; ++i)
            direction[i] =
                outsideScv.center()[i]
                - insideScv.center()[i];

        const Scalar centerDistance =
            direction.two_norm();

        if (!std::isfinite(centerDistance)
            || centerDistance <= Scalar(1.0e-20))
            return tensor;

        direction /= centerDistance;

        const auto& normal =
            scvf.unitOuterNormal();

        Scalar normalProjection = 0.0;

        for (int i = 0; i < dimWorld; ++i)
            normalProjection +=
                direction[i]*normal[i];

        Vector velocity(0.0);

        if (std::abs(normalProjection) > Scalar(1.0e-8))
        {
            const Scalar scale =
                qNormal/normalProjection;

            for (int i = 0; i < dimWorld; ++i)
                velocity[i] =
                    scale*direction[i];
        }
        else
        {
            // Safe fallback on an extremely non-orthogonal face.
            for (int i = 0; i < dimWorld; ++i)
                velocity[i] =
                    qNormal*normal[i];
        }

        const Scalar speed =
            velocity.two_norm();

        if (!std::isfinite(speed)
            || speed <= Scalar(1.0e-20))
            return tensor;

        static const Scalar alphaL =
            getParam<Scalar>(
                "Dispersion.LongitudinalDispersivity",
                0.50
            );

        static const Scalar alphaT =
            getParam<Scalar>(
                "Dispersion.TransverseDispersivity",
                0.05
            );

        // Scheidegger tensor
        //
        // D = alpha_T |q| I
        //   + (alpha_L-alpha_T) q q^T / |q|
        //
        // alphaL and alphaT remain exactly the HYDRUS values.
        for (int i = 0; i < dimWorld; ++i)
        {
            tensor[i][i] =
                alphaT*speed;

            for (int j = 0; j < dimWorld; ++j)
                tensor[i][j] +=
                    (alphaL-alphaT)
                    *velocity[i]
                    *velocity[j]
                    /speed;
        }

        return tensor;
    }

    static Scalar transmissibility_(
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf,
        const Matrix& tensor)
    {
        const auto insideIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideIdx);
        const auto& insideVolVars = elemVolVars[insideIdx];
        const Scalar ti = computeTpfaTransmissibility(
            fvGeometry, scvf, insideScv, tensor,
            insideVolVars.extrusionFactor()
        );

        if (scvf.boundary())
            return Extrusion::area(fvGeometry, scvf)*ti;

        const auto outsideIdx = scvf.outsideScvIdx();
        const auto& outsideScv = fvGeometry.scv(outsideIdx);
        const auto& outsideVolVars = elemVolVars[outsideIdx];
        const Scalar tj = -computeTpfaTransmissibility(
            fvGeometry, scvf, outsideScv, tensor,
            outsideVolVars.extrusionFactor()
        );

        // Mechanical dispersion must remain diffusive. On a
        // non-K-orthogonal triangular grid an anisotropic tensor can produce
        // a negative one-sided TPFA transmissibility. Reject such faces
        // instead of allowing an anti-diffusive component flux.
        if (!std::isfinite(ti)
            || !std::isfinite(tj)
            || ti <= Scalar(0.0)
            || tj <= Scalar(0.0))
            return Scalar(0.0);

        const Scalar denominator = ti + tj;
        if (!std::isfinite(denominator)
            || denominator <= Scalar(1.0e-30))
            return Scalar(0.0);

        const Scalar transmissibility =
            Extrusion::area(fvGeometry, scvf)
            *(ti*tj)/denominator;

        return std::isfinite(transmissibility)
               && transmissibility > Scalar(0.0)
               ? transmissibility
               : Scalar(0.0);
    }
};

} // namespace Dumux

#endif
