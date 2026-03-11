/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProbeFormBridge.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Form.hpp"
#include "MediumSystem.hpp"
#include "Probe.hpp"
#include "Snapshot.hpp"
#include "SpatialGrid.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

ProbeFormBridge::ProbeFormBridge(const Probe* probe, const Form* form)
{
    // remember the probe and form being bridged
    _probe = probe;
    _form = form;

    // find the simulation's spatial grid, if present
    if (probe->find<Configuration>()->hasMedium()) _grid = probe->find<MediumSystem>()->grid();

    // find the simulation's units system
    _units = probe->find<Units>();
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::writeQuantity(string fileid, string unit, string description, string projectedDescription,
                                    const Array& axis, string axisUnit, AddColumnDefinitions addColumnDefinitions,
                                    CompoundValueInCell valueInCell, WeightInCell weightInCell)
{
    _type = Type::GridCompoundAveraged;

    _fileid = fileid;
    _projectedFileid = _fileid;
    _unit = unit;
    _projectedUnit = _unit;
    _unitFactor = 1.;
    _projectedUnitFactor = 1.;
    _description = description;
    _projectedDescription = projectedDescription;
    _axis = axis;
    _axisUnit = axisUnit;
    _numValues = axis.size();

    _addColumnDefinitions = addColumnDefinitions;
    _compoundValueInCell = valueInCell;
    _weightInCell = weightInCell;

    _form->writeQuantity(this);
}

////////////////////////////////////////////////////////////////////

const SimulationItem* ProbeFormBridge::probe() const
{
    return _probe;
}

////////////////////////////////////////////////////////////////////

const SpatialGrid* ProbeFormBridge::grid() const
{
    return _grid;
}

////////////////////////////////////////////////////////////////////

const Units* ProbeFormBridge::units() const
{
    return _units;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // if the specified iteration index is positive, return a 3-digit decimal representation with leading zeros,
    // otherwise return the empty string
    string iterString(int iter)
    {
        string result;
        if (iter > 0)
        {
            result = std::to_string(iter);
            if (result.length() < 3) result.insert(0, 3 - result.length(), '0');
        }
        return result;
    }
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::prefix() const
{
    string result = _probe->itemName() + iterString(_probe->iter());
    if (!_fileid.empty()) result += "_" + _fileid;
    return result;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::projectedPrefix() const
{
    string result = _probe->itemName() + iterString(_probe->iter());
    if (!_projectedFileid.empty()) result += "_" + _projectedFileid;
    return result;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::unit() const
{
    return _unit;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::projectedUnit() const
{
    return _projectedUnit;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::description() const
{
    return _description;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::projectedDescription() const
{
    return _projectedDescription;
}

////////////////////////////////////////////////////////////////////

Array ProbeFormBridge::axis() const
{
    return _axis;
}

////////////////////////////////////////////////////////////////////

string ProbeFormBridge::axisUnit() const
{
    return _axisUnit;
}

////////////////////////////////////////////////////////////////////

int ProbeFormBridge::numValues() const
{
    return _numValues;
}

////////////////////////////////////////////////////////////////////

bool ProbeFormBridge::isVector() const
{
    return _type == Type::GridVectorAveraged || _type == Type::InputVector;
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::addColumnDefinitions(TextOutFile& outfile) const
{
    switch (_type)
    {
        case Type::GridScalarAccumulated:
        case Type::GridScalarAveraged:
        case Type::InputScalar:
        {
            outfile.addColumn(_description, _unit);
            break;
        }
        case Type::GridVectorAveraged:
        case Type::InputVector:
        {
            outfile.addColumn(_description + " x", _unit);
            outfile.addColumn(_description + " y", _unit);
            outfile.addColumn(_description + " z", _unit);
            break;
        }
        case Type::GridCompoundAccumulated:
        case Type::GridCompoundAveraged:
        case Type::InputCompound:
        {
            _addColumnDefinitions(outfile);
        }
    }
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::valuesInCell(int m, Array& values) const
{
    switch (_type)
    {
        case Type::GridScalarAccumulated:
        case Type::GridScalarAveraged:
        {
            values[0] = _scalarValueInCell(m) * _unitFactor;
            break;
        }
        case Type::GridVectorAveraged:
        {
            Vec v = _vectorValueInCell(m) * _unitFactor;
            values[0] = v.x();
            values[1] = v.y();
            values[2] = v.z();
            break;
        }
        case Type::GridCompoundAccumulated:
        {
            values = _compoundValueInCell(m);
            values *= _unitFactor;
            break;
        }
        case Type::GridCompoundAveraged:
        {
            values = _compoundValueInCell(m);
            break;
        }
        case Type::InputScalar:
        case Type::InputVector:
        case Type::InputCompound:
        {
            throw FATALERROR("Cannot retrieve values in given spatial cell from input model probe");
        }
    }
}

////////////////////////////////////////////////////////////////////
