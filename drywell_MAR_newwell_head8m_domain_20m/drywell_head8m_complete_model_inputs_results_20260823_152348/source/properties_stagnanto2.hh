// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef DUMUX_FLOODMAR_STAGNANT_O2_PROPERTIES_HH
#define DUMUX_FLOODMAR_STAGNANT_O2_PROPERTIES_HH

#include "properties.hh"
#include "floodmarstagnanto2localresidual.hh"

namespace Dumux::Properties {

namespace TTag {

struct FloodMarStagnantO2
{
    using InheritsFrom = std::tuple<FloodMar>;
};

} // namespace TTag

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FloodMarStagnantO2>
{
    using type = FloodMarStagnantO2LocalResidual<TypeTag>;
};

} // namespace Dumux::Properties

#endif
