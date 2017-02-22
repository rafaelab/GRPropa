#include "grpropa/module/Synchrotron.h"
#include "grpropa/Units.h"
#include "grpropa/Cosmology.h"
#include "grpropa/magneticField/MagneticField.h"

#include <limits>

namespace grpropa {

Synchrotron::Synchrotron(ref_ptr<MagneticField> field, double limit) 
{
    this->Bfield = field;
    this->limit = limit;
}

void Synchrotron::process(Candidate *c) const {
    // check if electrons / positrons
    if (std::abs(c->current.getId()) != 11)
        return;

    double z = c->getRedshift();
    double E = c->current.getEnergy() * (1 + z);
    double lf = c->current.getLorentzFactor() * (1 + z);
    Vector3d pos = c->current.getPosition();
    Vector3d b = Bfield->getField(pos);
    Vector3d v = c->current.getVelocity();
    double beta = c->current.getSpeed() / c_light;
    double step = c->getCurrentStep() / (1 + z);
    double B = b.getR() * pow(1 + z, 2);

    double vperp = v.cross(b).getR() / B;
    double RL = lf * mass_electron * vperp / (eplus * B);
    double dEdx = (2. / 3.) * pow(eplus, 2) * pow(beta * lf, 4) / pow(RL, 2); 

    double dE = std::abs(dEdx * step);
    dE = std::min(E, dE);
    double Enew = (E - dE) / (1 + z);

    if (Enew * (1 + z) > E)
        return;

    c->current.setEnergy(Enew);
    c->limitNextStep(limit * E / dEdx);
}


} // namespace grpropa
