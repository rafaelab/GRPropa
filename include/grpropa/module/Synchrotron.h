#ifndef GRPROPA_SYNCHROTRON_H
#define GRPROPA_SYNCHROTRON_H

#include "grpropa/Module.h"
#include "grpropa/magneticField/MagneticField.h"

namespace grpropa {


/**
 @class Synchrotron
 @brief Energy losses for electrons due to synchrotron emission.

 This implementation follows the one of the Elmag code [Kachelriess et al. 10.1016/j.cpc.2011.12.025].\n
 The synchrotron emission is considered a continuous energy loss process. \n
 Losses are calculated step by step using dE/dX from the parametrization from :\n
   V.N. Baier, V.M. Katkov and V.M. Strakhovenko, “Electromagnetic processes at high energies in oriented single crystals”, World Scientific (1998).\n

 */
class Synchrotron: public Module {
private:
    ref_ptr<MagneticField> Bfield;
    double CriticalB;
public:
    Synchrotron(ref_ptr<MagneticField> field, double Bcr = 4.14e9);
	void process(Candidate *candidate) const;
};


} // namespace grpropa

#endif // GRPROPA_REDSHIFT_H
