#ifndef MPDLAMBDAPOL_H
#define MPDLAMBDAPOL_H

/// \ingroup physics
/// \class MpdLambda
/// \brief Lambda reconstructed object
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - 22-03-2023
/// Modified to include extra information for polarization analysis (Elizaveta Nazarova)

#include "TObject.h"

class MpdLambdaPol {

 public:
   MpdLambdaPol() {;}
   MpdLambdaPol(Float_t imassh, Float_t ipth, Float_t iph, Float_t ietah, Float_t iyh, Float_t ichi2h, Float_t idisth, 
	     Float_t ipath, Float_t iangle, Float_t *ietas, Float_t *imcthetas, Float_t *ithetas, Float_t *imcphis,
       Float_t *iphis, Float_t *imcps, Float_t *ips, Float_t *ipts, Float_t *ichi2s, Float_t *idcas,
	     Float_t idca, Float_t ic2pv, Float_t iomega1, Float_t iomega2, Float_t icosA, Float_t icosAmc, 
       Float_t ipolarhx, Float_t ipolarhy, Float_t ipolarhz, Float_t iphi_star, Float_t iphi_star_MC, Float_t iphi_Lam, 
       Float_t *ic2s, Int_t *iorigs, Int_t *iqs, Int_t *ilayMx, Int_t impdg);
   virtual ~MpdLambdaPol() {}

 public:
   Float_t massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas[2], ps[2], mcps[2], pts[2];
   Float_t chi2s[2], c2s[2], dcas[2], mcthetas[2], thetas[2], mcphis[2], phis[2];
   Float_t dca, c2pv, omega1, omega2, cosA, cosAmc, polarhx, polarhy, polarhz, phi_star, phi_star_MC, phi_Lam;
   Int_t origs[2], qs[2], layMx[2], mpdg;

   // private:

   ClassDef(MpdLambdaPol, 1);
};
#endif
