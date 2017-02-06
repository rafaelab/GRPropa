#ifndef GRPROPA_PAIRPRODUCTION_H
#define GRPROPA_PAIRPRODUCTION_H


#include "grpropa/Module.h"
#include "grpropa/PhotonBackground.h"

namespace grpropa {

/**
 @class PairProduction
 @brief Pair production of photons with background photons.

 This implementation follows the one of the Elmag code [Kachelriess et al. 10.1016/j.cpc.2011.12.025].\n
 This module simulates electron-pair production as an stochastic process.\n
 Several photon fields can be selected, although the dominant one is the infrared.\n
 By default, the module limits the step size to 10% of the energy loss length of the particle.
 */
class PairProduction: public Module {
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
    
public:
    PairProduction(PhotonField photonField = CMB, double thinning = 0., double limit = 0.1, double nMaxIterations = 1000);

    void setPhotonField(PhotonField photonField);
    void setLimit(double limit);
    void setThinning(double thinning);
    void setMaxNumberOfIterations(double nMaxIterations);
    void initTableBackgroundEnergy(std::string filename);
    void initRate(std::string filename);
    void process(Candidate *candidate) const;
    double centerOfMassEnergy2(double E, double e, double mu) const; 
    double energyFraction(double E, double z) const;
    double lossLength(int id, double en, double z) const;
    void performInteraction(Candidate *candidate) const;
};

} // namespace grpropa


#endif // GRPROPA_PAIRPRODUCTION_H
