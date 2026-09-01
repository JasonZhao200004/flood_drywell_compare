// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef DUMUX_FLOODMAR_CURRENT_PHASEGUARD_V2_FACE_DISP_NO_GAS_O2_ADV_HH
#define DUMUX_FLOODMAR_CURRENT_PHASEGUARD_V2_FACE_DISP_NO_GAS_O2_ADV_HH

#include "properties_oldboundary_phaseguard_v2_facedisp.hh"
#include "floodmar_noGasO2Adv_current_localresidual.hh"

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
        FloodMarCurrentNoGasO2AdvLocalResidual<TypeTag>;
};

} // namespace Dumux::Properties

#endif
