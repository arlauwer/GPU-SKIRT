/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConvergenceInfoProbe.hpp"
#include "Configuration.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // writes two output lines with input and gridded values, respectively
    void writeValues(TextOutFile& out, double t, double g, double factor, string unit)
    {
        string tail;
        if (!unit.empty()) tail += " " + unit;
        out.writeLine("  Input:   " + StringUtils::toString(t * factor, 'g', 6) + tail);
        if (t > 0) tail += " (" + StringUtils::toString(100. * (g - t) / t, 'f', 2) + " %)";
        out.writeLine("  Gridded: " + StringUtils::toString(g * factor, 'g', 6) + tail);
    }

    // outputs the convergence info for the given material type
    void writeConvergenceForMaterialType(TextOutFile& out, MediumSystem* ms, MaterialMix::MaterialType type,
                                         string name)
    {
        Units* units = ms->find<Units>();
        int numMedia = ms->numMedia();
        int numCells = ms->numCells();

        // output header
        out.writeLine("");
        out.writeLine("------ " + name + " ------");

        // calculate the true and gridded total mass
        double mass_t = 0.;
        double mass_g = 0.;
        for (int h = 0; h != numMedia; ++h)
        {
            if (ms->isMaterialType(type, h))
            {
                mass_t += ms->media()[h]->mass();
                for (int m = 0; m != numCells; ++m) mass_g += ms->massDensity(m, h) * ms->volume(m);
            }
        }

        // output the total mass
        out.writeLine("");
        out.writeLine("Total mass");
        writeValues(out, mass_t, mass_g, units->omass(1.), units->umass());
    }

    // outputs statistics on optical depth
    void writeOpticalDepthStatistics(TextOutFile& out, MediumSystem* ms, double lambda)
    {
        Units* units = ms->find<Units>();

        // calculate diagonal optical depth for all spatial cells
        int numCells = ms->numCells();
        Array tauV(numCells);
        for (int m = 0; m != numCells; ++m) tauV[m] = ms->grid()->diagonal(m) * ms->opacityExt(lambda, m);

        // calculate statistics on optical depth
        double tauavg = tauV.sum() / numCells;
        double taumin = tauV.min();
        double taumax = tauV.max();
        const int numBins = 500;
        vector<int> countV(numBins + 1);
        for (double tau : tauV)
        {
            int index = max(0, min(numBins, static_cast<int>((tau - taumin) / (taumax - taumin) * numBins)));
            countV[index]++;
        }
        int count = 0;
        int index = 0;
        for (; index < numBins; index++)
        {
            count += countV[index];
            if (count > 0.9 * numCells) break;
        }
        double tau90 = taumin + index * (taumax - taumin) / numBins;

        // write the statistics on optical depth to the file
        out.writeLine("");
        out.writeLine("------ Cell statistics ------");
        out.writeLine("");
        out.writeLine("Optical depth of cell diagonal at " + StringUtils::toString(units->owavelength(lambda), 'g')
                      + " " + units->uwavelength());
        out.writeLine("  Largest: " + StringUtils::toString(taumax, 'g'));
        out.writeLine("  Average: " + StringUtils::toString(tauavg, 'g'));
        out.writeLine("  90% <  : " + StringUtils::toString(tau90, 'g'));
    }
}

////////////////////////////////////////////////////////////////////

void ConvergenceInfoProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // create a text file
        TextOutFile out(this, itemName() + "_convergence", "spatial grid convergence info");
        out.writeLine("Spatial grid convergence information");
        out.writeLine("------------------------------------");

        // write info for each material type, if present
        using Type = MaterialMix::MaterialType;
        if (ms->hasDust()) writeConvergenceForMaterialType(out, ms, Type::Dust, "Dust");
        if (ms->hasElectrons()) writeConvergenceForMaterialType(out, ms, Type::Electrons, "Electrons");
        if (ms->hasGas()) writeConvergenceForMaterialType(out, ms, Type::Gas, "Gas");

        // write statistics on optical depth
        writeOpticalDepthStatistics(out, ms, wavelength());

        // write footer
        out.writeLine("");
        out.writeLine("------------------------------------");
    }
}

////////////////////////////////////////////////////////////////////
