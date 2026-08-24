// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef DUMUX_FLOODMAR_HYDRUSMATCH_PROPERTIES_HH
#define DUMUX_FLOODMAR_HYDRUSMATCH_PROPERTIES_HH
#include "properties.hh"
#include "floodmar_hydrusmatch_stagnanto2_localresidual.hh"
namespace Dumux::Properties {
namespace TTag {
struct FloodMarHydrusMatch { using InheritsFrom = std::tuple<FloodMar>; };
}
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FloodMarHydrusMatch>
{
    using type = FloodMarHydrusMatchStagnantO2LocalResidual<TypeTag>;
};
}
#endif
