#include "grpropa/module/InverseCompton.h"
#include "grpropa/Random.h"
#include "grpropa/Units.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace grpropa {

InverseCompton::InverseCompton(PhotonField photonField, double limit, double ethr) {
	setPhotonField(photonField);
	this->limit = limit;
    this->Ethr = ethr;
}

void InverseCompton::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("InverseCompton: CMB");
		initRate(getDataPath("MFP_IC.dat"));
		initTableBackgroundEnergy(getDataPath("CMBcum.dat"));
		break;
	case EBL:  // default: Kneiske '04 IRB model
	case EBL_Kneiske10:
		setDescription("InverseCompton: IRB Kneiske '10 (lower limit)");
		initRate(getDataPath("epp_IRB_Kneiske10.txt"));
		break;
	case EBL_Franceschini08:
		setDescription("InverseCompton: IRB Franceschini '08");
		initRate(getDataPath("epp_IRB_Franceschini08.txt"));
		break;
	default:
		throw std::runtime_error("InverseCompton: unknown photon background");
	}
}


void InverseCompton::setLimit(double limit) {
	this->limit = limit;
}

void InverseCompton::setThresholdEnergy(double Ethr) {
    this->Ethr = Ethr;
}

void InverseCompton::initRate(std::string filename) {
	
    std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("InverseCompton: could not open file " + filename);

	// clear previously loaded interaction rates
	tabEnergy.clear();
	tabRate.clear();
  
    // For now approximating ICS interaction rate as constant.
    // Should be replaced by the correct values.
    
    /*int n = 50;
	for (int i=0; i<n; i++){
		double en = pow(10, 10+.1*i);
		double t = 1.2e11;
		tabEnergy.push_back(en * eV);
		tabRate.push_back(1/(c_light*t));
	}*/
    
    
	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy.push_back(a * eV);
				tabRate.push_back(1 / (b * Mpc));
                //std::cout << b << "  "  << a << std::endl;
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
    
}

void InverseCompton::initTableBackgroundEnergy(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("InverseCompton: could not open file " + filename);
	tabPhotonEnergy.clear();
	tabProb.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile){
				tabPhotonEnergy.push_back(pow(10,a) * eV);
				tabProb.push_back(b);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

double InverseCompton::energyFraction(double E, double z) const {
	/* 
		Returns the fraction of energy of the incoming electron taken by the
		outgoing electron.
		See ELMAG code for this implementation.
	    Kachelriess et al., Comp. Phys. Comm 183 (2012) 1036

       Note: E is the correct value and should not be multiplied by (1+z).
             The redshift is used in this function only to draw the energy of the 
             background photon.
	*/
	Random &random = Random::instance();

	// drawing energy of background photon according to number density (integral)
	double e = interpolate(random.rand(), tabProb, tabPhotonEnergy);
    e *= (1 + z);
    double ethr = ethr * (1 + z);

	// kinematics
	double mu = random.randUniform(-1, 1);
	double s = centerOfMassEnergy2(E, e, mu);
	double ymin = pow(mass_electron * c_squared, 2) / s;
    double eps = ethr / E;
	double ymax = 1 - eps; 
	double y;
	double r = random.rand();
	double gb = 0;
	do {
		y = ymin * pow(ymax / ymin,r);
		double f1 = (1 + y*y) / 2;
		double f2 = 2 * ymin * (y - ymin) * (1 - y) / (y * pow(1 - ymin, 2));
		double gb = (f1 - f2);
		if (random.rand() < gb)
			break;
	} while(r < gb);
	//std::cout << s / pow(eV,2) << " " << ymin << " " << 1-y << " " << e/eV << " " << E/eV << std::endl;
	return y;
}

double InverseCompton::energyLossBelowThreshold(double E, double z, double step) const {
	Random &random = Random::instance();

	// drawing energy of background photon according to number density (integral)
	double e = interpolate(random.rand(), tabProb, tabPhotonEnergy);
    e *= (1 + z);
    double ethr = Ethr * (1 + z);

    const double ThomsonCS = 6.65e24;

    // kinematics
	double mu = random.randUniform(-1, 1);
	double s = centerOfMassEnergy2(E, e, mu);
	double ymin = pow(mass_electron * c_squared, 2) / s;
    double eps = ethr / E;
	double ymax = 1 - eps; // should be 1 - eps
    double a0 = (1 - ymax) / (1 - ymin);
    double a1 = (log(1 / ymax) / (1 - ymax) - 1) * (1 - 4 * ymin * (1 + 2 * ymin) / pow(1 - ymin, 2));
    double a2 = (1 - ymax) * (1 + 2 * ymax) / 6;
    double a3 = 2 * ymin * (1 + 2 * ymin / ymax) * (1 - ymax) / pow(1 - ymin, 2);
    double dE = 0.75 * ThomsonCS * ymin * a0 * (a1 + a2 + a3);
    return dE * step;
}

double InverseCompton::centerOfMassEnergy2(double E, double e, double mu) const {
	double beta = sqrt(1 - pow(mass_electron*c_squared, 2) / (E * E));
	return pow(mass_electron*c_squared, 2) + 2 * E * e * (1 - beta * mu);
}


double InverseCompton::lossLength(int id, double en, double z) const {

	en *= (1 + z);
	if (en < tabEnergy.front())
		return std::numeric_limits<double>::max(); // below energy threshold

	double rate;
	if (en < tabEnergy.back())
		rate = interpolate(en, tabEnergy, tabRate); // interpolation
	else
		rate = tabRate.back() * pow(en / tabEnergy.back(), -0.6); // extrapolation

	rate *= pow(1 + z, 3);
	return 1. / rate;
}

void InverseCompton::process(Candidate *c) const {
	// execute the loop at least once for limiting the next step
	double step = c->getCurrentStep();
	int id = c->current.getId();
    double E = c->current.getEnergy();
	double z = c->getRedshift();
    E *= (1 + z);

    // only electrons / positrons allowed
	if (fabs(id) != 11)
		return; 

    //if (E > Ethr) {
    do {
  	    //double rate = interpolate(E, tabEnergy, tabRate);
        double rate = 1/ lossLength(id, E / (1 + z), z);
	    // cosmological scaling, rate per comoving distance)
        rate *= pow(1 + z, 2);

		Random &random = Random::instance();
	    double randDistance = -log(random.rand()) / rate;

		// check if an interaction occurs in this step
		if (step < randDistance) {
		   	// limit next step to a fraction of the mean free path
		  	c->limitNextStep(limit / rate);
		    return;
		}

		performInteraction(c);

		// repeat with remaining steps
		step -= randDistance;
	} while (step > 0);
    //} 
    /*else { // continuous loss approximation for soft photons
        double E = c->current.getEnergy();
        double dE = energyLossBelowThreshold(E, z, step);
        c->current.setEnergy(E - dE);
    }*/
}

void InverseCompton::performInteraction(Candidate *candidate) const {
	
	double en = candidate->current.getEnergy();
	double z = candidate->getRedshift();
	double y = energyFraction(en, z);
	candidate->current.setEnergy(en * y);
    candidate->setActive(false);
	candidate->addSecondary(22, en * (1 - y));
	//std::cout << "22\t" << en/eV << "\t" << y << std::endl;

}


} // namespace grpropa
