#ifndef GRPROPA_OUTPUTTXT_H
#define GRPROPA_OUTPUTTXT_H

#include "grpropa/Module.h"

#include <bitset>
#include <stdexcept>

namespace grpropa {

/**
 @class Output
 @brief Configurable output base class.
 */
class Output: public Module {
protected:
    double lengthScale, energyScale;
    std::bitset<64> fields;
    bool oneDimensional;
    mutable size_t count;
    
    void modify();

public:
    enum OutputColumn {
        WeightColumn,
        CosmicTimeColumn,
        TrajectoryLengthColumn,
        RedshiftColumn,
        CurrentIdColumn,
        CurrentEnergyColumn,
        CurrentPositionColumn,
        CurrentDirectionColumn,
        SourceIdColumn,
        SourceEnergyColumn,
        SourcePositionColumn,
        SourceDirectionColumn,
        CreatedIdColumn,
        CreatedEnergyColumn,
        CreatedPositionColumn,
        CreatedDirectionColumn
    };

    enum OutputType {
        Trajectory1D,
        Trajectory3D,
        Event1D,
        Event3D,
        Everything
    };

    Output();
    Output(OutputType outputtype);
    
    void setEnergyScale(double scale);
    void setLengthScale(double scale);

    void setOutputType(OutputType outputtype);
    void set(OutputColumn field, bool value);
    void enable(OutputColumn field);
    void disable(OutputColumn field);
    void enableAll();
    void disableAll();
    void set1D(bool value);
    size_t getCount() const;

    void process(Candidate *) const;
};



} // namespace grpropa

#endif // GRPROPA_OUTPUTTXT_H
