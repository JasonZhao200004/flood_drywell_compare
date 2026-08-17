// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef FLOODMAR_TRANSIENT_CCTPFA_DISPERSION_FLUX_HH
#define FLOODMAR_TRANSIENT_CCTPFA_DISPERSION_FLUX_HH

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
            elemFluxVarsCache, phaseIdx
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
        const ElementFluxVariablesCache& elemFluxVarsCache,
        const int phaseIdx)
    {
        static_assert(dimWorld == 2,
                      "The project-local velocity reconstruction is 2-D.");

        // Fit v from n_f dot v = q_f/(phi*S), using every face of this cell.
        Matrix normalMatrix(0.0);
        Vector rhs(0.0);

        for (const auto& face : scvfs(fvGeometry))
        {
            // A Neumann/Robin face has no physical outside volume variables.
            // Its prescribed flux remains owned by problem.hh, so reconstruct
            // the cell velocity from internal faces only.
            if (face.boundary())
                continue;

            const Scalar area = Extrusion::area(fvGeometry, face);
            if (area < 1.0e-30)
                continue;

            const auto& vi = elemVolVars[face.insideScvIdx()];
            const auto& vo = elemVolVars[face.outsideScvIdx()];
            const Scalar mobility = 0.5*(vi.mobility(phaseIdx)
                                       + vo.mobility(phaseIdx));
            const Scalar porosity = 0.5*(vi.porosity() + vo.porosity());
            const Scalar saturation = std::max(
                Scalar(1.0e-12),
                Scalar(0.5)*(vi.saturation(phaseIdx) + vo.saturation(phaseIdx))
            );

            // AdvectionType::flux is the integrated Darcy driving flux before
            // multiplication by mobility. qNormal is Darcy velocity [m/s].
            const Scalar qNormal =
                AdvectionType::flux(problem, element, fvGeometry, elemVolVars,
                                    face, phaseIdx, elemFluxVarsCache)
                * mobility / area;
            const Scalar poreNormal = qNormal
                / std::max(Scalar(1.0e-12), porosity*saturation);

            const auto& n = face.unitOuterNormal();
            for (int i = 0; i < dimWorld; ++i)
            {
                rhs[i] += n[i]*poreNormal;
                for (int j = 0; j < dimWorld; ++j)
                    normalMatrix[i][j] += n[i]*n[j];
            }
        }

        Vector velocity(0.0);
        const Scalar determinant = normalMatrix[0][0]*normalMatrix[1][1]
                                 - normalMatrix[0][1]*normalMatrix[1][0];
        if (std::abs(determinant) > 1.0e-30)
        {
            velocity[0] = ( normalMatrix[1][1]*rhs[0]
                          - normalMatrix[0][1]*rhs[1])/determinant;
            velocity[1] = (-normalMatrix[1][0]*rhs[0]
                          + normalMatrix[0][0]*rhs[1])/determinant;
        }

        static const Scalar alphaL =
            getParam<Scalar>("Dispersion.LongitudinalDispersivity", 0.50);
        static const Scalar alphaT =
            getParam<Scalar>("Dispersion.TransverseDispersivity", 0.05);

        const Scalar speed = velocity.two_norm();
        Matrix tensor(0.0);
        if (speed < 1.0e-20)
            return tensor;

        for (int i = 0; i < dimWorld; ++i)
        {
            tensor[i][i] = alphaT*speed;
            for (int j = 0; j < dimWorld; ++j)
                tensor[i][j] += (alphaL-alphaT)
                               *velocity[i]*velocity[j]/speed;
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

        if (ti*tj <= 0.0)
            return 0.0;
        return Extrusion::area(fvGeometry, scvf)*(ti*tj)/(ti+tj);
    }
};

} // namespace Dumux

#endif
