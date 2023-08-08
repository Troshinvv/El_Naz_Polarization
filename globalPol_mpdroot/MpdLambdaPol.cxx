/// \ingroup physics
/// \class MpdLambda
/// \brief Lambda reconstructed object
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - 22-03-2023
/// Modified to include extra information for polarization analysis (Elizaveta Nazarova)

#include "MpdLambdaPol.h"

ClassImp(MpdLambdaPol);

//-----------------------------------------------------------------------------------------------

MpdLambdaPol::MpdLambdaPol(Float_t imassh, Float_t ipth, Float_t iph, Float_t ietah, Float_t iyh, Float_t ichi2h,
                           Float_t idisth, Float_t ipath, Float_t iangle, Float_t *ietas, Float_t *imcthetas,
                           Float_t *ithetas, Float_t *imcphis, Float_t *iphis, Float_t *imcps, Float_t *ips,
                           Float_t *ipts, Float_t *ichi2s, Float_t *idcas, Float_t idca, Float_t ic2pv, Float_t iomega1,
                           Float_t iomega2, Float_t icosA, Float_t icosAmc, Float_t ipolarhx, Float_t ipolarhy,
                           Float_t ipolarhz, Float_t iphi_star, Float_t iphi_star_MC, Float_t iphi_Lam, Float_t *ic2s,
                           Int_t *iorigs, Int_t *iqs, Int_t *ilayMx, Int_t impdg)
   : massh(imassh), pth(ipth), ph(iph), etah(ietah), yh(iyh), chi2h(ichi2h), disth(idisth), path(ipath), angle(iangle),
     dca(idca), c2pv(ic2pv), omega1(iomega1), omega2(iomega2), cosA(icosA), cosAmc(icosAmc), polarhx(ipolarhx),
     polarhy(ipolarhy), polarhz(ipolarhz), phi_star(iphi_star), phi_star_MC(iphi_star_MC), phi_Lam(iphi_Lam),
     mpdg(impdg)
{

   for (Int_t j = 0; j < 2; ++j) {
      etas[j]     = ietas[j];
      ps[j]       = ips[j];
      pts[j]      = ipts[j];
      chi2s[j]    = ichi2s[j];
      dcas[j]     = idcas[j];
      c2s[j]      = ic2s[j];
      origs[j]    = iorigs[j];
      qs[j]       = iqs[j];
      layMx[j]    = ilayMx[j];
      phis[j]     = iphis[j];
      thetas[j]   = ithetas[j];
      mcphis[j]   = imcphis[j];
      mcthetas[j] = imcthetas[j];
      mcps[j]     = imcps[j];
   }
}

//-----------------------------------------------------------------------------------------------