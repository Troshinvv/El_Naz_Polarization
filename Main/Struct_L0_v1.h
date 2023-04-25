#ifndef STRUCT_L0_V1_H
#define STRUCT_L0_V1_H
#include <stdlib.h>
#include <vector>

//structure for L0
struct L0 {
public:
  L0() {}
  L0(Float_t imassh, Float_t ipth, Float_t iph, Float_t ietah, Float_t iyh, Float_t ichi2h, Float_t idisth, 
     Float_t ipath, Float_t iangle, Float_t *ietas, Float_t *imcthetas, Float_t *ithetas, Float_t *imcphis,
     Float_t *iphis, Float_t *imcps, Float_t *ips, Float_t *ipts, Float_t *ichi2s, Float_t *idcas, 
     Float_t idca, Float_t ic2pv, Float_t iomega1, Float_t iomega2, Float_t icosA, Float_t icosAmc, Float_t ipolarhx, Float_t ipolarhy, Float_t ipolarhz, Float_t iphi_star, Float_t iphi_star_MC, Float_t iphi_Lam,
     Int_t *iorigs, Int_t *iqs, Int_t *ilayMx, Int_t ievNo) : 
    massh(imassh), pth(ipth), ph(iph), etah(ietah), yh(iyh), chi2h(ichi2h), disth(idisth), path(ipath),
    angle(iangle), dca(idca), c2pv(ic2pv), omega1(iomega1), omega2(iomega2), cosA(icosA), cosAmc(icosAmc), polarhx(ipolarhx), polarhy(ipolarhy), polarhz(ipolarhz), phi_star(iphi_star), phi_star_MC(iphi_star_MC), phi_Lam(iphi_Lam), evNo(ievNo) {
    for (Int_t j = 0; j < 2; ++j) {
      etas[j] = ietas[j];
      ps[j] = ips[j];
      pts[j] = ipts[j];
      chi2s[j] = ichi2s[j];
      dcas[j] = idcas[j];
      origs[j] = iorigs[j];
      qs[j] = iqs[j];
      layMx[j] = ilayMx[j];
      phis[j] = iphis[j];
      thetas[j] = ithetas[j];
      mcphis[j] = imcphis[j];
      mcthetas[j] = imcthetas[j];
      mcps[j] = imcps[j];
  }}
  Float_t massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas[2], ps[2], mcps[2], pts[2], chi2s[2], dcas[2],
    mcthetas[2], thetas[2], mcphis[2], phis[2], cosA, cosAmc, polarhx, polarhy, polarhz, phi_star, phi_star_MC, phi_Lam;
  Float_t dca, c2pv, omega1, omega2;
  Int_t origs[2], qs[2], layMx[2], evNo;
};
//structure for Xi
struct Xi {
public:
  Xi() {}
  Xi(Float_t imassh, Float_t ipth, Float_t iph, Float_t ietah, Float_t iyh, Float_t ichi2h, Float_t idisth, 
     Float_t ipath, Float_t iangle, Float_t ic2pv, Float_t idca, Float_t *ietas, Float_t *ips, 
     Float_t *ipts, Float_t *ichi2s, Float_t *idcas, Float_t imassL, Float_t ichi2L, Float_t idistL,
     Float_t ipathL,  Float_t iangL,  Float_t *ichi2sL, Float_t *idcasL, Float_t imassL1, Float_t iomega1, Float_t iomega2,
     Float_t *iomegaL, Int_t *iorigs, Int_t *iqs, Int_t *ilayMx, Int_t ievNo) : 
    massh(imassh), pth(ipth), ph(iph), etah(ietah), yh(iyh), chi2h(ichi2h), disth(idisth), path(ipath),
    angle(iangle), c2pv(ic2pv), dca(idca), massL(imassL), chi2L(ichi2L), distL(idistL), pathL(ipathL), 
    angL(iangL), massL1(imassL1), omega1(iomega1), omega2(iomega2), evNo(ievNo) {
    for (Int_t j = 0; j < 2; ++j) {
      etas[j] = ietas[j];
      ps[j] = ips[j];
      pts[j] = ipts[j];
      chi2s[j] = ichi2s[j];
      chi2sL[j] = ichi2sL[j];
      dcasL[j] = idcasL[j];
      dcas[j] = idcas[j];
      origs[j] = iorigs[j];
      qs[j] = iqs[j];
      layMx[j] = ilayMx[j];
    }
    for (Int_t j = 0; j < 3; ++j) omegaL[j] = iomegaL[j];
  }
  Float_t massh, pth, ph, etah, yh, chi2h, disth, path, angle, c2pv, dca, etas[2], ps[2], pts[2], chi2s[2], dcas[2];
  Float_t massL, chi2L, distL, pathL, angL, chi2sL[2], dcasL[2], massL1, omega1, omega2, omegaL[3];
  Int_t origs[2], qs[2], layMx[2], evNo;
};
#endif
#ifdef __ROOTCLING__
#pragma link C++ class L0+;
#pragma link C++ class std::vector<L0>+;
#pragma link C++ class Xi+;
#pragma link C++ class std::vector<Xi>+;
#pragma link C++ class std::vector<pair<float,float> >+;
#pragma link C++ class std::tuple<float,float,float>+;
#pragma link C++ class std::vector<tuple<float,float,float> >+;
#endif
