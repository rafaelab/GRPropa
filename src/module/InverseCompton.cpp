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
        redshiftDependence = false;
        setDescription("Inverse Compton: CMB");
        initRate(getDataPath("ICS-CMB.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-CMB.txt"));
        break;
    case EBL:  // default: Gilmore '12 IRB model
    case EBL_Gilmore12:
        redshiftDependence = true;
        setDescription("Inverse Compton: EBL Gilmore et al. 2012");
        initRate(getDataPath("ICS-EBL_Gilmore12.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Gilmore12.txt"));
        break;
    case EBL_Dominguez11:
        redshiftDependence = true;
        setDescription("Inverse Compton: EBL Dominguez et al. 2011");
        initRate(getDataPath("ICS-EBL_Dominguez11.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Dominguez11.txt"));
        break;
    case EBL_Dominguez11_UL:
        redshiftDependence = true;
        setDescription("Inverse Compton: EBL Dominguez et al. 2011 (upper limit)");
        initRate(getDataPath("ICS-EBL_Dominguez11_upper.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Dominguez11_upper.txt"));
        break;
    case EBL_Dominguez11_LL:
        redshiftDependence = true;
        setDescription("Inverse Compton: EBL Dominguez et al. 2011 (lower limit)");
        initRate(getDataPath("ICS-EBL_Dominguez11_lower.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Dominguez11_lower.txt"));
        break;
    case EBL_Finke10:
        redshiftDependence = true;
        setDescription("Inverse Compton: EBL Finke et al. 2010");
        initRate(getDataPath("ICs-EBL_Finke10.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Finke10.txt"));
        break;
    case EBL_Kneiske10:
        redshiftDependence = true;
        setDescription("Inverse Compton: EBL Kneiske & Dole 2010 (lower limit)");
        initRate(getDataPath("ICS-EBL_Kneiske10.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Kneiske10.txt"));
        break;
    case EBL_Franceschini08:
        redshiftDependence = true;
        setDescription("Inverse Compton: EBL Franceschini et al. 2008");
        initRate(getDataPath("ICS-EBL_Franceschini08.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-EBL_Franceschini08.txt"));
        break;
    case CRB:
    case CRB_Protheroe96:
        redshiftDependence = false;
        setDescription("Inverse Compton: CRB Protheroe & Biermann 1996");
        initRate(getDataPath("ICS-CRB_Protheroe96.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-CRB_Protheroe96.txt"));
        break;
    case CRB_ARCADE2:
        redshiftDependence = false;
        setDescription("Inverse Compton: CRB ARCADE2 2010");
        initRate(getDataPath("ICS-CRB_ARCADE2.txt"));
        initTableBackgroundEnergy(getDataPath("photonProbabilities-CRB_ARCADE2.txt"));
        break;
    default:
        throw std::runtime_error("Inverse Compton: unknown photon background");
    }
}


void InverseCompton::setLimit(double limit) {
    this->limit = limit;
}

void InverseCompton::setThresholdEnergy(double Ethr) {
    this->Ethr = Ethr;
}

void InverseCompton::initRate(std::string filename) {
    
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
            nc = 33;
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
            double redshifts[] = {0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.9};
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

    // for (int i=0; i<tabRate.size(); i++) std::cout << tabRate[i] * Mpc << std::endl;
    // for (int i=0; i<tabEnergy.size(); i++) std::cout << tabEnergy[i] /eV << std::endl;
    
}

void InverseCompton::initTableBackgroundEnergy(std::string filename) {
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
                infile >> b >> a;
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
        int nl = 400; // number of lines (probabilities)

        if (photonField == EBL_Finke10) nc = 33;
        else if (photonField == EBL_Gilmore12) nc = 20;
        else if (photonField == EBL_Kneiske10) nc = 5;
        else if (photonField == EBL_Franceschini08) nc = 11;
        else if (photonField == EBL_Dominguez11 || photonField == EBL_Dominguez11_UL || photonField == EBL_Dominguez11_LL) nc = 18;
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

        for (int j=0; j<nl; j++) 
            tabProb.push_back(entries[0][j]);
        for (int i=1; i<nc+1; i++){
            for (int j=0; j<nl; j++){
                tabPhotonEnergy.push_back(entries[i][j] * eV);
            }
        }
        infile.close();

        // test
        // for (int i=0; i<tabPhotonEnergy.size(); i++) std::cout << tabPhotonEnergy[i] / eV << std::endl;
        // for (int i=0; i<tabProb.size(); i++) std::cout << tabProb[i] << std::endl;
        // for (int i=0; i<tabRedshift.size(); i++) std::cout << tabRedshift[i] << std::endl;
    } // conditional: redshift dependent
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
        y = ymin * pow(ymax / ymin, r);
        double f1 = (1 + y * y) / 2;
        double f2 = 2 * ymin * (y - ymin) * (1 - y) / (y * pow(1 - ymin, 2));
        double gb = (f1 - f2);
        if (random.rand() < gb)
            break;
    } while(r < gb);
    //std::cout << s / pow(eV,2) << " " << ymin << " " << 1-y << " " << e/eV << " " << E/eV << std::endl;
    if (y > 0 && y < 1)
        return y;
    else
        return -1;
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
    if (y > 0 && y < 1) {
        candidate->current.setEnergy(en * y);
        candidate->setActive(false);
        candidate->addSecondary(22, en * (1 - y));
    }
}


} // namespace grpropa
