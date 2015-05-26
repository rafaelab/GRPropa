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
        initRate(getDataPath("PP-CMB.txt"));
        initTableBackgroundEnergy(getDataPath("cumDensities-CMB.txt"));
        break;
    case EBL:  // default: Gilmore '12 IRB model
    case EBL_Gilmore12:
        redshiftDependence = true;
        setDescription("PairProduction: Gilmore et al. 2012");
        initRate(getDataPath("PP-EBL_Gilmore12.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Gilmore12.txt"));
        break;
    case EBL_Dominguez11:
        redshiftDependence = true;
        setDescription("PairProduction: Dominguez et al. 2011");
        initRate(getDataPath("PP-EBL_Dominguez11.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Dominguez11.txt"));
        break;
    case EBL_Dominguez11_UL:
        redshiftDependence = true;
        setDescription("PairProduction: Dominguez et al. 2011 (upper limit)");
        initRate(getDataPath("PP-EBL_Dominguez11_upper.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Dominguez11_upper.txt"));
        break;
    case EBL_Dominguez11_LL:
        redshiftDependence = true;
        setDescription("PairProduction: Dominguez et al. 2011 (lower limit)");
        initRate(getDataPath("PP-EBL_Dominguez11_lower.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Dominguez11_lower.txt"));
        break;
    case EBL_Finke10:
        redshiftDependence = true;
        setDescription("PairProduction: EBL Finke et al. 2010");
        initRate(getDataPath("PP-EBL_Finke10.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Finke10.txt"));
        break;
    case EBL_Kneiske10:
        redshiftDependence = true;
        setDescription("PairProduction: EBL Kneiske & Dole 2010 (lower limit)");
        initRate(getDataPath("PP-EBL_Kneiske10.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Kneiske10.txt"));
        break;
    case EBL_Franceschini08:
        redshiftDependence = true;
        setDescription("PairProduction: EBL Franceschini et al. 2008");
        initRate(getDataPath("PP-EBL_Franceschini08.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Franceschini08.txt"));
        break;
    case CRB:
    case CRB_Protheroe96:
        redshiftDependence = false;
        setDescription("PairProduction: CRB Protheroe & Biermann 1996");
        initRate(getDataPath("PP-CRB_Protheroe96.txt"));
        initTableBackgroundEnergy(getDataPath("cumDensities-CRB_Protheroe96.txt"));
        break;
    case CRB_ARCADE2:
        redshiftDependence = false;
        setDescription("PairProduction: CRB ARCADE2 2010");
        initRate(getDataPath("PP-CRB_ARCADE2.txt"));
        initTableBackgroundEnergy(getDataPath("cumDensities-CRB_ARCADE2.txt"));
        break;
    default:
        throw std::runtime_error("PairProduction: unknown photon background");
    }
}

void PairProduction::setLimit(double limit) {
    this->limit = limit;
}

void PairProduction::initRate(std::string filename) {
    
    if (redshiftDependence == false) {
        std::ifstream infile(filename.c_str());
        if (!infile.good())
            throw std::runtime_error("PairProduction: could not open file " + filename);
   
        // clear previously loaded interaction rates
        tabEnergy.clear();
        tabRate.clear();

        while (infile.good()) {
            if (infile.peek() != '#') {
                double a, b;
                infile >> a >> b;
                if (infile) {
                    tabEnergy.push_back(a * eV);
                    tabRate.push_back(b / Mpc);
                }
    		}
    		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
        }
        infile.close();
    } else { // Rates for EBL
        std::ifstream infile(filename.c_str());
 	    if (!infile.good())
		    throw std::runtime_error("PairProduction: could not open file " + filename);
  
	    // clear previously loaded interaction rates
    	tabEnergy.clear();
    	tabRate.clear();
        tabRedshift.clear();

        // size of vector is predefined and depends on the model
        int nc; // number of columns (redshifts + one column for energy)
        int nl = 81; // number of lines (energies)
        std::vector<double> redshifts;
        if (photonField == EBL_Finke10) {
            nc = 23;
            double redshifts[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.09, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 4.99};
            for (int k=0; k<nc; k++) 
                tabRedshift.push_back(redshifts[k]);
        }
        else if (photonField == EBL_Gilmore12) {
            nc = 20;
            double redshifts[] = {0, 0.015, 0.025, 0.044, 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0};
            for (int k=0; k<nc; k++) 
                tabRedshift.push_back(redshifts[k]);
        }
        else if (photonField == EBL_Dominguez11 || photonField == EBL_Dominguez11_UL || photonField == EBL_Dominguez11_LL) {
            nc = 20;
            double redshifts[] = {0, 0.015, 0.025, 0.044, 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0};
            for (int k=0; k<nc; k++) 
                tabRedshift.push_back(redshifts[k]);
        }
        else if (photonField == EBL_Kneiske10) {
            nc = 5;
            double redshifts[] = {0.0, 0.1, 0.3, 0.8, 2.0};
            for (int k=0; k<nc; k++) 
                tabRedshift.push_back(redshifts[k]);
     	}
        else if (photonField == EBL_Franceschini08) {
        	nc = 11;
        	double redshifts[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
        	for (int k=0; k<nc; k++) 
        		tabRedshift.push_back(redshifts[k]);
        }
        else { 
        	throw std::runtime_error("EBL model not defined for redshift dependent treatment (or not defined at all).");
        }

        double entries[nc+1][nl];
        int j = 0;
        while (!infile.eof()) {
        	for (int i=0; i<nc+1; i++) {
        		double entry = 0;
        		infile >> entry;
        		entries[i][j] = entry;
        	}
        	j++;
        }

        for (int j=0; j<nl; j++) 
        	tabEnergy.push_back(entries[0][j] * eV);
        for (int i=1; i<nc+1; i++){
        	for (int j=0; j<nl; j++){
        		tabRate.push_back(entries[i][j] / Mpc);
        	}
        }
        infile.close();

    } // conditional: redshift dependent
}

void PairProduction::initTableBackgroundEnergy(std::string filename) {

 
	if (redshiftDependence == false) {
	    std::ifstream infile(filename.c_str());
	    if (!infile.good())
		    throw std::runtime_error("PairProduction: could not open file " + filename);
   
	    // clear previously loaded interaction rates
    	tabPhotonEnergy.clear();
	   	tabProb.clear();

    	while (infile.good()) {
    		if (infile.peek() != '#') {
    			double a, b;
    			infile >> a >> b;
    			if (infile) {
    				tabPhotonEnergy.push_back(a * eV);
    				tabProb.push_back(b);
    			}
    		}
    		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    	}
    	infile.close();
    } else { // Rates for EBL
        std::ifstream infile(filename.c_str());
 	    if (!infile.good())
		    throw std::runtime_error("PairProduction: could not open file " + filename);
  
	    // clear previously loaded interaction rates
    	tabPhotonEnergy.clear();
    	tabProb.clear();

        // size of vector is predefined and depends on the model
        int nc; // number of columns (redshifts + one column for energy)
        int nl = 250; // number of lines (energies)

		if (photonField == EBL_Finke10) nc = 500;
       	else if (photonField == EBL_Gilmore12) nc = 20;
        else if (photonField == EBL_Kneiske10) nc = 5;
        else if (photonField == EBL_Franceschini08)	nc = 11;
        else throw std::runtime_error("EBL model not defined for redshift dependent treatment (or not defined at all).");

        double entries[nc+1][nl];
        int j = 0;
        while (!infile.eof()) {
        	for (int i=0; i<nc+1; i++) {
        		double entry = 0;
        		infile >> entry;
        		entries[i][j] = entry;
        	}
        	j++;
        }
        for (int i=0; i<nc-1; i++){
        	for (int j=0; j<)

        }

        // for (int j=0; j<nl; j++) 
        // 	tabPhotonEnergy.push_back(entries[0][j] * eV);
        // for (int i=1; i<nc+1; i++){
        // 	for (int j=0; j<nl; j++){
        // 		tabProb.push_back(entries[i][j]);
        // 	}
        // }
        infile.close();

    } // conditional: redshift dependent
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
	//double e = interpolate(random.rand(), tabProb, tabPhotonEnergy);
	double rnd = random.rand();
	// double rnd = 1e-4;
	double e = 0;
	if (redshiftDependence == true) {
		e = interpolate2d(z, rnd, tabRedshift, tabProb, tabPhotonEnergy);
		std::cout << z << "\t" << rnd << "\t" << e / eV << std::endl;
	 } else {
		e = (1 + z) * interpolate(rnd, tabProb, tabPhotonEnergy);
		// std::cout << z << "\t" << rnd << "\t" << e / eV << std::endl;
	}

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
		} while(r < gb);
		if (random.rand() > 0.5) 
			y = 1 - y;
        if (y <= 0 || y >= 1) {
        	y = 0.5;
        	std::cout << "WARNING: energy fraction of produced pair out of range: 0<y<1" << std::endl;  
        }
        // std::cout << z << " " << E << " " << e << " " << s << " " << y << " " << std::endl;
		return y;
	} else {
		return 0;
	}

     
}

double PairProduction::centerOfMassEnergy2(double E, double e, double mu) const {
	return 2 * E * e * (1 - mu);
}


double PairProduction::lossLength(int id, double en, double z) const {
	
	if (id != 22)
		return std::numeric_limits<double>::max(); // no pair production on uncharged particles

	if (en < tabEnergy.front())
		return std::numeric_limits<double>::max(); // below energy threshold

	double rate;
    if (redshiftDependence == false) {
		en *= (1 + z);
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
		}
	}

	return 1. / rate;
}


void PairProduction::process(Candidate *c) const {

	// execute the loop at least once for limiting the next step
	double step = c->getCurrentStep();
	do {
		int id = c->current.getId();
		if (id != 22) 
			return; // only photons allowed

		double en = c->current.getEnergy();
		double z = c->getRedshift();
        double rate = 1 / lossLength(id, en, z);

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
        // This occurs to a very small fraction of events, so this workaround 
        // doesn't affect significantly the results.
	    candidate->addSecondary(11, en * .5);
        candidate->addSecondary(-11, en * .5);
    }
}

} // namespace grpropa
