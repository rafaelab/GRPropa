#ifndef GRPROPA_PHOTONBACKGROUND_H
#define GRPROPA_PHOTONBACKGROUND_H

namespace grpropa {

// Photon fields
// The default IRB model is that of Finke et al. 2010
enum PhotonField {
  CMB, 
  EBL, EBL_Finke10, EBL_Kneiske10, EBL_Franceschini08, EBL_Gilmore12, EBL_Dominguez11, EBL_Dominguez11_UL, EBL_Dominguez11_LL, 
  CRB, CRB_Protheroe96, CRB_ARCADE2
};

// Returns list of available photon background models
void listOfPhotonBackgroundModels();


} // namespace grpropa

#endif // GRPROPA_PHOTONBACKGROUND_H
