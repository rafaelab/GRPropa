#ifndef GRPROPA_PHOTONBACKGROUND_H
#define GRPROPA_PHOTONBACKGROUND_H

namespace grpropa {

// Photon fields
// The default IRB model is that of Finke et al. 2010
enum PhotonField {
  CMB, IRB, IRB_Finke10, IRB_Kneiske04, IRB_Kneiske10, IRB_Franceschini08, IRB_Stecker05, IRB_withRedshift_Kneiske04
};

// Returns overall comoving scaling factor
double photonFieldScaling(PhotonField photonField, double z);

} // namespace grpropa

#endif // GRPROPA_PHOTONBACKGROUND_H
