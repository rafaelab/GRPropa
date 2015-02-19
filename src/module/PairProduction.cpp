#include "grpropa/module/PairProduction.h"
#include "grpropa/Random.h"
#include "grpropa/Units.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace grpropa {

PairProduction::PairProduction(PhotonField photonField, double limit) {
	setPhotonField(photonField);
	this->limit = limit;
}

void PairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
        redshiftDependence = false;
		setDescription("PairProduction: CMB");
		initRate(getDataPath("MFP-PP-CMB.dat"));
		initTableBackgroundEnergy(getDataPath("CMBcum.dat"));
		break;
	case IRB:  // default: Finke '10 IRB model
    case IRB_Finke10:
        redshiftDependence = true;
		setDescription("PairProduction: IRB Finke '10");
		initRate(getDataPath("MFP-PP-CIB-z.dat"));
		initTableBackgroundEnergy(getDataPath("CIBCumModel1.dat"));
		break;
    case IRB_Kneiske04:
        redshiftDependence = false;
		setDescription("PairProduction: IRB test");
		initRate(getDataPath("MFP-PP-CIB-z0.00.dat"));
		initTableBackgroundEnergy(getDataPath("CIBCumModel1.dat"));
		break;
	default:
		throw std::runtime_error("PairProduction: unknown photon background");
	}
}

void PairProduction::setLimit(double limit) {
	this->limit = limit;
}

void PairProduction::initRate(std::string filename) {
    
    // Rates for CMB   
    if (photonField == CMB) { 
	    std::ifstream infile(filename.c_str());
	    if (!infile.good())
		    throw std::runtime_error("PairProduction: could not open file " + filename);
   
	    // clear previously loaded interaction rates
    	tabEnergy.clear();
        tabLogEnergy.clear();
    	tabRate.clear();

    	while (infile.good()) {
    		if (infile.peek() != '#') {
    			double a, b;
    			infile >> a >> b;
    			if (infile) {
    				tabEnergy.push_back(a * eV);
                    tabLogEnergy.push_back(log10(a * eV));
    				tabRate.push_back(b / Mpc);
    			}
    		}
    		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    	}
    	infile.close();
    } else {
        if (redshiftDependence == false) {
    	    std::ifstream infile(filename.c_str());
    	    if (!infile.good())
    		    throw std::runtime_error("PairProduction: could not open file " + filename);
     
    	    // clear previously loaded interaction rates
        	tabEnergy.clear();
            tabLogEnergy.clear();
        	tabRate.clear();

        	while (infile.good()) {
        		if (infile.peek() != '#') {
        			double a, b;
        			infile >> a >> b;
        			if (infile) {
        				tabEnergy.push_back(a * eV);
                        tabLogEnergy.push_back(log10(a * eV));
        				tabRate.push_back(b / Mpc);
        			}
        		}
        		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
        	}
        	infile.close();
        } else { // Redshift dependent case
    	    std::ifstream infile(filename.c_str());
     	    if (!infile.good())
    		    throw std::runtime_error("PairProduction: could not open file " + filename);
      
    	    // clear previously loaded interaction rates
        	tabEnergy.clear();
        	tabRate.clear();
            tabLogEnergy.clear();


            // size of vector is predefined.
            int nc = 95; // number of columns
            int nl = 501; // number of lines
            std::vector< std::vector<double> > rates(nl, std::vector<double>(nc));
            //int j = 0; // current line
            int i = 0;
            while (!infile.eof()){
                if (i == 0) {
                    double dummy, e;
                    infile >> dummy;
                    for (int j=0; j<nc; j++) {
                        infile >> e;
                        tabEnergy.push_back(e * eV);
                        tabLogEnergy.push_back(log10(e * eV));
                    }
                } else {
                    double z;
                    infile >> z;
                    tabRedshift.push_back(z);
                     for (int j=0; j<nc; j++) {
                        double r;
                        infile >> r;
                        tabRate.push_back(r / Mpc);
                    }
                }
                i++;
            }
        	infile.close();
        }
    }

}

void PairProduction::initTableBackgroundEnergy(std::string filename) {

	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("PairProduction: could not open file " + filename);
	tabPhotonEnergy.clear();
	tabProb.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile){
				tabPhotonEnergy.push_back(a * eV);
				tabProb.push_back(b);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

double PairProduction::energyFraction(double E, double z) const {
	/* 
		Returns the fraction of energy of the incoming photon taken by the
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

	// kinematics
	double mu = random.randUniform(-1, 1);	
	double s = 2 * E * e * (1 - mu);
	if (4 * pow(mass_electron * c_squared, 2) / s < 1) {
		double beta = sqrt(1 - 4 * pow(mass_electron * c_squared, 2) / s);
		double ymin = (1 - beta) / 2;
		double y = 0;
		double r = random.rand();
		double gb = 0;
		do {
			y = 0.5 * pow(2 * ymin, random.rand());
			double pf = 1 / (1 + 2 * beta * beta * (1 - beta * beta));
			double f1 = y * y / (1 - y);
			double f2 = 1 - y + (1 - beta * beta) / (1 - y);
			double f3 = pow(1 - beta * beta, 2) / (4. * y * pow(1 - y, 2) );
			double gb = pf * (f1 + f2 - f3);
            //std::cout << "s, beta, ymin, gb, pf, f1, f2, f3 " << s << " " << beta << " " << ymin << " " << gb << " " << pf " " << f1 << " " << f2 << " " << f3 << std::endl;
            //std::cout << gb << "  "  << pf << "  " << f1 << "  "  << f2 << "  "  << f3 << std::endl;
		} while(r < gb);
		if (random.rand() > 0.5) 
			y = 1 - y;
        if (y<=0 || y>=1) std::cout << "warning: energy fraction of pair out of range: 0<y<1" << std::endl;  
            //std::cout << "s, beta, ymin, gb " << s << " " << beta << " " << ymin << " " << gb << std::endl;
		return y;
	}
	else return 0;


     
}

double PairProduction::centerOfMassEnergy2(double E, double e, double mu) const {
	return 2 * E * e * (1 - mu);
}



double PairProduction::lossLength(int id, double en, double z) const {
	
	if (id != 22)
		return std::numeric_limits<double>::max(); // no pair production on uncharged particles

	en *= (1 + z);
	if (en < tabEnergy.front())
		return std::numeric_limits<double>::max(); // below energy threshold

	double rate;
    if (redshiftDependence == false) {
	    if (en < tabEnergy.back())
		    rate = interpolate(en, tabEnergy, tabRate); // interpolation
	    else
		    rate = tabRate.back() * pow(en / tabEnergy.back(), -0.6); // extrapolation
    	rate *= pow(1 + z, 3) * photonFieldScaling(photonField, z);
    } else {
        if (en < tabEnergy.back()) {
		    rate = interpolate2d(z, en, tabRedshift, tabEnergy, tabRate); // interpolation
	    } else {
		    rate = tabRate.back() * pow(en / tabEnergy.back(), -0.6); // extrapolation
    	    //rate *= pow(1 + z, 3) * photonFieldScaling(photonField, z);
        }    
    }
	return 1. / rate;
}


void PairProduction::process(Candidate *c) const {

	// execute the loop at least once for limiting the next step
	double step = c->getCurrentStep();
	do {
		int id = c->current.getId();
		if (id != 22) {
			return; // only photons allowed
		}
		double en = c->current.getEnergy();
		double z = c->getRedshift();
        double rate = 0;
        rate = 1 / lossLength(id, en, z);


		Random &random = Random::instance();
		double randDistance = -log(random.rand()) / rate;
        //std::cout << step/kpc << "\t" << randDistance/Mpc << std::endl;

		// check if an interaction occurs in this step
		if (step < randDistance) {
			// limit next step to a fraction of the mean free path
			c->limitNextStep(limit / rate);
    		return;
		}

        //std::cout << step/kpc << "\t" << randDistance/kpc << "\t" << rate*Mpc << "\t" << 1 / (rate*Mpc) << "\t" << en/eV <<  "\t"  << limit/rate/Mpc << std::endl;

		performInteraction(c);

		// repeat with remaining steps
		step -= randDistance;
	} while (step > 0);
}

void PairProduction::performInteraction(Candidate *candidate) const {
	
	double en = candidate->current.getEnergy();
	double z = candidate->getRedshift();
	double y = energyFraction(en, z);
    candidate->setActive(false);
    if (y > 0 && y < 1) {
	    candidate->addSecondary(11, en * y);
        candidate->addSecondary(-11, en * (1 - y));
    } else { 
        // Fixes bug of y=0. Problem should be understood.
        // This occurs to a very small fraction of events, 
        // so this workaround is an excellent approximation.
	    candidate->addSecondary(11, en * .5);
        candidate->addSecondary(-11, en * .5);
    }
}

} // namespace grpropa
