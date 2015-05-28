#ifndef GRPROPA_REDSHIFT_H
#define GRPROPA_REDSHIFT_H

#include "grpropa/Module.h"

namespace grpropa {

/**
 @class Redshift
 @brief Updates redshift and applies adiabatic energy loss according to the travelled distance.
 */
class Redshift: public Module {
public:
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class FutureRedshift
 @brief Updates redshift and applies adiabatic energy loss according to the travelled distance.
 */
class FutureRedshift: public Module {
public:
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};

} // namespace grpropa

#endif // GRPROPA_REDSHIFT_H
