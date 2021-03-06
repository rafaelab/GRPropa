#ifndef GRPROPA_INVERSECOMPTON_H
#define GRPROPA_INVERSECOMPTON_H

#include "grpropa/Module.h"
#include "grpropa/Units.h" 
#include "grpropa/PhotonBackground.h"

#include <vector>

namespace grpropa {

/**
 @class InverseCompton
 @brief Inverse Compton scattering of electrons off background photons.

 This module simulates inverse Compton scattering as continuous energy loss below Ethr (10 GeV by default).\n
 This implementation follows the one of the Elmag code [Kachelriess et al. 10.1016/j.cpc.2011.12.025] \n
 Several photon fields can be selected, although CMB is the dominant one.\n
 For now supports only electrons/positrons, but corresponding effect for muons may be included in the future,
 in spite of the fact that it is virtually negligible.\n
 By default, the module limits the step size to 10% of the energy loss length of the particle.
 */
class InverseCompton: public Module {
private:
    PhotonField photonField;

    std::vector<double> tabEnergy; /* tabulated energy [eV] */
    std::vector<double> tabRate; /* tabulated rate [1/Mpc] */
    std::vector<double> tabRedshift; /* tabulated redshifts for z dependence of the IRB */
    std::vector<double> tabPhotonEnergy; /* background photon energy*/
    std::vector<double> tabProb; /* cumulative probability for background photon. */

    double thinning; /* number of secondaries to be tracked; if 1 only one secondary is tracked */
    double limit; /* fraction of energy loss length to limit the next step */
    double nMaxIterations; /* maximum number of attempts to sample s in energy fraction */
    bool redshiftDependence; /* whether EBL model is redshift-dependent */
    double Ethr;  /* energy loss due to the emission of soft photons for E<Ethr */

public:
    InverseCompton(PhotonField photonField = CMB, double thinning = 0, double limit = 0.1, double Ethr = 1e6 * eV, double nMaxInteractions = 1000);

    void setPhotonField(PhotonField photonField);
    void setLimit(double limit);
    void setThinning(double thinning);
    void setThresholdEnergy(double Ethr);
    void setMaxNumberOfIterations(double nMaxIterations);
    void initRate(std::string filename);
    void initTableBackgroundEnergy(std::string filename);
    void process(Candidate *candidate) const;
    double lossLength(int id, double lf, double z) const;
    double energyLossBelowThreshold(double E, double z, double step) const; 
    double centerOfMassEnergy2(double E, double e, double mu) const; 
    double energyFraction(double E, double z) const;
    void performInteraction(Candidate *candidate) const;
};

} // namespace grpropa

#endif // GRPROPA_INVERSECOMPTON_H
