#ifndef GRPROPA_SYNCHROTRON_H
#define GRPROPA_SYNCHROTRON_H

#include "grpropa/Module.h"
#include "grpropa/magneticField/MagneticField.h"

namespace grpropa {


/**
 @class Synchrotron
 @brief Energy losses for electrons due to synchrotron emission.

 See "Classical electrodynamics" by J. D. Jackson, eq. 14.31, pg. 667.
 This implementation only accounts for the radiated power. 
 The spectrum of radiated particles will be implemented in later versions of the code.
 */
class Synchrotron: public Module {
private:
    ref_ptr<MagneticField> Bfield;
    double limit;
public:
    Synchrotron(ref_ptr<MagneticField> field, double limit = 0.1);
    void process(Candidate *candidate) const;
};


} // namespace grpropa

#endif // GRPROPA_SYNCHROTRON_H
