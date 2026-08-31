// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef DUMUX_FLOODMAR_HYDRUSLIKE_NOGASO2ADV_PROPERTIES_HH
#define DUMUX_FLOODMAR_HYDRUSLIKE_NOGASO2ADV_PROPERTIES_HH

#include "properties_oldboundary_phaseguard_v2_facedisp.hh"
#include "floodmar_hydruslike_noGasO2Adv_localresidual.hh"

namespace Dumux::Properties {

namespace TTag {

struct FloodMarOldBoundaryPhaseGuardV2FaceDispNoGasO2Adv
{
    using InheritsFrom =
        std::tuple<FloodMarOldBoundaryPhaseGuardV2FaceDisp>;
};

} // namespace TTag

template<class TypeTag>
struct LocalResidual<
    TypeTag,
    TTag::FloodMarOldBoundaryPhaseGuardV2FaceDispNoGasO2Adv
>
{
    using type =
        FloodMarHydrusLikeNoGasO2AdvLocalResidual<TypeTag>;
};

} // namespace Dumux::Properties

#endif
