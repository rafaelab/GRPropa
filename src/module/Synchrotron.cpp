#include "grpropa/module/Synchrotron.h"
#include "grpropa/Units.h"
#include "grpropa/Cosmology.h"
#include "grpropa/magneticField/MagneticField.h"

#include <limits>

namespace grpropa {

Synchrotron::Synchrotron(ref_ptr<MagneticField> field, double Bcr) 
{
    this->Bfield = field;
    this->CriticalB = Bcr;
}

void Synchrotron::process(Candidate *c) const {
	// check if electrons / positrons
	if (abs(c->current.getId()) != 11 )
		return;

	double z = c->getRedshift();
	double E = c->current.getEnergy();
    double gamma = c->current.getLorentzFactor();
    Vector3d pos = c->current.getPosition();
    Vector3d p = c->current.getMomentum();    
    Vector3d b = Bfield->getField(pos);
    double B = b.getR();

    //E *= 1 + z;
    /*double pi = 3.1415926;
    double hbar = h_planck / (2 * pi);
    double chi = (p.cross(b)).getR() / ( mass_electron * mass_electron * c_light * CriticalB );
    double dE = pow(mass_electron * c_squared, 2) / pow( 1 + 4.8 * (1 + chi) * log(1 + 1.7 * chi)  + 3.44 * chi * chi, 2/3 ) * (eV / (hbar * c_light)) / 137.;*/ // from Elmag
    double dE = 1.058e-14 * pow(B * gamma, 2) * (1 - 1 / pow(gamma, 2)); // from Longarir
    dE = dE * c->getCurrentStep();
    dE = std::min(E, dE); // energy loss should not be larger than energy of particle
	c->current.setEnergy(E * (1 - dE));
}



} // namespace grpropa
