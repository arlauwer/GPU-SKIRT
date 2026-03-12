/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultCutsForm.hpp"
#include "MediumSystem.hpp"
#include "PlanarCutsForm.hpp"

////////////////////////////////////////////////////////////////////

DefaultCutsForm::DefaultCutsForm(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void DefaultCutsForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    // the size in pixels (in each spatial direction) for the default cuts
    constexpr int Np = 1024;

    // the dimension of the medium system and the extent of the spatial grid
    auto ms = find<MediumSystem>();
    Box box = ms->grid()->boundingBox();

    // output cuts depending on the dimension of the medium system
    PlanarCutsForm::writePlanarCut(bridge, 1, 1, 0, box, 0., 0., 0., Np, Np, Np);  // xy
    PlanarCutsForm::writePlanarCut(bridge, 1, 0, 1, box, 0., 0., 0., Np, Np, Np);  // xz
    PlanarCutsForm::writePlanarCut(bridge, 0, 1, 1, box, 0., 0., 0., Np, Np, Np);  // yz
}

////////////////////////////////////////////////////////////////////
