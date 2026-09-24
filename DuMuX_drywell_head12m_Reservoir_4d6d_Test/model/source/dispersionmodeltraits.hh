// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef FLOODMAR_DISPERSION_MODEL_TRAITS_HH
#define FLOODMAR_DISPERSION_MODEL_TRAITS_HH

namespace Dumux {

template<class BaseTraits>
struct FloodMarDispersionModelTraits : public BaseTraits
{
    static constexpr bool enableCompositionalDispersion() { return true; }
};

} // namespace Dumux

#endif
