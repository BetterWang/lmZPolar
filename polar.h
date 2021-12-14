#include "const.h"
#include "HiEvtPlaneListPbPb.h"
#include "HiEvtPlaneListpPb.h"
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <vector>
#include <array>
#include <complex>

using namespace std;

int NCent = 0;
const double *CentBins = nullptr;

inline double roundPi(double phi) {
    while ( phi > TMath::Pi() ) phi-=TMath::Pi();
    while ( phi <-TMath::Pi() ) phi+=TMath::Pi();
    return phi;
}

struct PolarEvent {
    // Declaration of leaf types
    vector<double>  *EPOrg      = 0;
    vector<double>  *EPOrgSin   = 0;
    vector<double>  *EPOrgCos   = 0;
    vector<double>  *EPFlat     = 0;
    vector<double>  *EPFlatSin  = 0;
    vector<double>  *EPFlatCos  = 0;
    vector<double>  *pdgId      = 0;
    vector<double>  *pTrkPt     = 0;
    vector<double>  *pTrkPx     = 0;
    vector<double>  *pTrkPy     = 0;
    vector<double>  *pTrkPz     = 0;
    vector<double>  *pTrkEta    = 0;
    vector<double>  *nTrkPt     = 0;
    vector<double>  *nTrkPx     = 0;
    vector<double>  *nTrkPy     = 0;
    vector<double>  *nTrkPz     = 0;
    vector<double>  *nTrkEta    = 0;
    vector<double>  *pPxCM      = 0;
    vector<double>  *pPyCM      = 0;
    vector<double>  *pPzCM      = 0;
    vector<double>  *nPxCM      = 0;
    vector<double>  *nPyCM      = 0;
    vector<double>  *nPzCM      = 0;
    vector<double>  *pt         = 0;
    vector<double>  *phi        = 0;
    vector<double>  *eta        = 0;
    vector<double>  *rapidity   = 0;
    vector<double>  *mass       = 0;
    Double_t        Cent        = 0;
    Double_t        vz          = 0;
    Double_t        RunId       = 0;
    Double_t        EventId     = 0;
    Double_t        Lumi        = 0;
};

struct PolarData {
    PolarData(std::string name, bool bFwd = false):
    hHF2pO_    (NCent),
    hHF2mO_    (NCent),
    hHF2pF_    (NCent),
    hHF2mF_    (NCent),
    hHF3pO_    (NCent),
    hHF3mO_    (NCent),
    hHF3pF_    (NCent),
    hHF3mF_    (NCent),
    hNLambda_  (NCent),
    hLambdaEta_(NCent),
    hLambdaMass_  (NCent, vector<TH1D*>(NpT)),
    hLamBarMass_  (NCent, vector<TH1D*>(NpT)),
    hLambdaP2Cos_ (NCent, vector<TH1D*>(NpT)),
    hLambdaP2Cos2_(NCent, vector<TH1D*>(NpT)),
    hLamBarP2Cos_ (NCent, vector<TH1D*>(NpT)),
    hLamBarP2Cos2_(NCent, vector<TH1D*>(NpT)),
    hLambdaM2Cos_ (NCent, vector<TH1D*>(NpT)),
    hLambdaM2Cos2_(NCent, vector<TH1D*>(NpT)),
    hLamBarM2Cos_ (NCent, vector<TH1D*>(NpT)),
    hLamBarM2Cos2_(NCent, vector<TH1D*>(NpT)),
    hLambdaP3Cos_ (NCent, vector<TH1D*>(NpT)),
    hLamBarP3Cos_ (NCent, vector<TH1D*>(NpT)),
    hLambdaM3Cos_ (NCent, vector<TH1D*>(NpT)),
    hLamBarM3Cos_ (NCent, vector<TH1D*>(NpT)),
    hLambdaF2MCos_(NCent, vector<TH1D*>(NpT)),
    hLambdaF2Cos2_(NCent, vector<TH1D*>(NpT)),
    hLambdaF3MCos_(NCent, vector<TH1D*>(NpT)),
    hLambdaF2PCos_(NCent, vector<TH1D*>(NpT)),
    hLambdaF3PCos_(NCent, vector<TH1D*>(NpT)),
    hLamBarF2MCos_(NCent, vector<TH1D*>(NpT)),
    hLamBarF2Cos2_(NCent, vector<TH1D*>(NpT)),
    hLamBarF3MCos_(NCent, vector<TH1D*>(NpT)),
    hLamBarF2PCos_(NCent, vector<TH1D*>(NpT)),
    hLamBarF3PCos_(NCent, vector<TH1D*>(NpT)),
    hLambdaPV2M_(NCent, vector<TH1D*>(NpT)),
    hLambdaMV2P_(NCent, vector<TH1D*>(NpT)),
    hLambdaFV2M_(NCent, vector<TH1D*>(NpT)),
    hLambdaFV2P_(NCent, vector<TH1D*>(NpT)),
    hLamBarPV2M_(NCent, vector<TH1D*>(NpT)),
    hLamBarMV2P_(NCent, vector<TH1D*>(NpT)),
    hLamBarFV2M_(NCent, vector<TH1D*>(NpT)),
    hLamBarFV2P_(NCent, vector<TH1D*>(NpT)),
    hLambdaPV3M_(NCent, vector<TH1D*>(NpT)),
    hLambdaMV3P_(NCent, vector<TH1D*>(NpT)),
    hLambdaFV3M_(NCent, vector<TH1D*>(NpT)),
    hLambdaFV3P_(NCent, vector<TH1D*>(NpT)),
    hLamBarPV3M_(NCent, vector<TH1D*>(NpT)),
    hLamBarMV3P_(NCent, vector<TH1D*>(NpT)),
    hLamBarFV3M_(NCent, vector<TH1D*>(NpT)),
    hLamBarFV3P_(NCent, vector<TH1D*>(NpT))
    {
        name_ = name;
        bForward_ = bFwd;
        if ( bpPb ) {
            hCent_ = new TH1D((name_+"hCent").c_str(), "hCent", 600, 0, 600);
        } else {
            hCent_ = new TH1D((name_+"hCent").c_str(), "hCent", 200, 0, 200);
        }

        hReHFp2HFm2_    = new TH1D((name_+"hReHFp2HFm2").c_str(),    "hReHFp2HFm2",    NCent, CentBins);
        hReHFp2TrkMid2_ = new TH1D((name_+"hReHFp2TrkMid2").c_str(), "hReHFp2TrkMid2", NCent, CentBins);
        hReHFm2TrkMid2_ = new TH1D((name_+"hReHFm2TrkMid2").c_str(), "hReHFm2TrkMid2", NCent, CentBins);
        hImHFp2HFm2_    = new TH1D((name_+"hImHFp2HFm2").c_str(),    "hImHFp2HFm2",    NCent, CentBins);
        hImHFp2TrkMid2_ = new TH1D((name_+"hImHFp2TrkMid2").c_str(), "hImHFp2TrkMid2", NCent, CentBins);
        hImHFm2TrkMid2_ = new TH1D((name_+"hImHFm2TrkMid2").c_str(), "hImHFm2TrkMid2", NCent, CentBins);

        hReHFp3HFm3_    = new TH1D((name_+"hReHFp3HFm3").c_str(),    "hReHFp3HFm3",    NCent, CentBins);
        hReHFp3TrkMid3_ = new TH1D((name_+"hReHFp3TrkMid3").c_str(), "hReHFp3TrkMid3", NCent, CentBins);
        hReHFm3TrkMid3_ = new TH1D((name_+"hReHFm3TrkMid3").c_str(), "hReHFm3TrkMid3", NCent, CentBins);
        hImHFp3HFm3_    = new TH1D((name_+"hImHFp3HFm3").c_str(),    "hImHFp3HFm3",    NCent, CentBins);
        hImHFp3TrkMid3_ = new TH1D((name_+"hImHFp3TrkMid3").c_str(), "hImHFp3TrkMid3", NCent, CentBins);
        hImHFm3TrkMid3_ = new TH1D((name_+"hImHFm3TrkMid3").c_str(), "hImHFm3TrkMid3", NCent, CentBins);

        for ( int cent = 0; cent < NCent; cent++ ) {
            hHF2pO_[cent] = new TH1D((name_+Form("hHF2pO_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2mO_[cent] = new TH1D((name_+Form("hHF2mO_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2pF_[cent] = new TH1D((name_+Form("hHF2pF_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2mF_[cent] = new TH1D((name_+Form("hHF2mF_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 72, -TMath::Pi()/2., TMath::Pi()/2.);

            hHF3pO_[cent] = new TH1D((name_+Form("hHF3pO_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 72, -TMath::Pi()/3., TMath::Pi()/3.);
            hHF3mO_[cent] = new TH1D((name_+Form("hHF3mO_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 72, -TMath::Pi()/3., TMath::Pi()/3.);
            hHF3pF_[cent] = new TH1D((name_+Form("hHF3pF_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 72, -TMath::Pi()/3., TMath::Pi()/3.);
            hHF3mF_[cent] = new TH1D((name_+Form("hHF3mF_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 72, -TMath::Pi()/3., TMath::Pi()/3.);

            hNLambda_[cent] = new TH1D((name_+Form("hNLambda_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 20, 0, 20);
            hLambdaEta_[cent] = new TH1D((name_+Form("hLambdaEta_%i", cent)).c_str(), Form("CentBin %i-%i", int(CentBins[cent]), int(CentBins[cent+1])), 50, -2.5, 2.5);

            for ( int ipt = 0; ipt < NpT; ipt++ ) {
                hLambdaMass_[cent][ipt] = new TH1D((name_+Form("hLambdaMass_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarMass_[cent][ipt] = new TH1D((name_+Form("hLamBarMass_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaP2Cos_ [cent][ipt] = new TH1D((name_+Form("hLambdaP2Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaP2Cos2_[cent][ipt] = new TH1D((name_+Form("hLambdaP2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarP2Cos_ [cent][ipt] = new TH1D((name_+Form("hLamBarP2Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarP2Cos2_[cent][ipt] = new TH1D((name_+Form("hLamBarP2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaM2Cos_ [cent][ipt] = new TH1D((name_+Form("hLambdaM2Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaM2Cos2_[cent][ipt] = new TH1D((name_+Form("hLambdaM2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarM2Cos_ [cent][ipt] = new TH1D((name_+Form("hLamBarM2Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarM2Cos2_[cent][ipt] = new TH1D((name_+Form("hLamBarM2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaP3Cos_ [cent][ipt] = new TH1D((name_+Form("hLambdaP3Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarP3Cos_ [cent][ipt] = new TH1D((name_+Form("hLamBarP3Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaM3Cos_ [cent][ipt] = new TH1D((name_+Form("hLambdaM3Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarM3Cos_ [cent][ipt] = new TH1D((name_+Form("hLamBarM3Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaF2MCos_[cent][ipt] = new TH1D((name_+Form("hLambdaF2MCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaF2Cos2_[cent][ipt] = new TH1D((name_+Form("hLambdaF2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaF3MCos_[cent][ipt] = new TH1D((name_+Form("hLambdaF3MCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaF2PCos_[cent][ipt] = new TH1D((name_+Form("hLambdaF2PCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaF3PCos_[cent][ipt] = new TH1D((name_+Form("hLambdaF3PCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLamBarF2MCos_[cent][ipt] = new TH1D((name_+Form("hLamBarF2MCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarF2Cos2_[cent][ipt] = new TH1D((name_+Form("hLamBarF2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarF3MCos_[cent][ipt] = new TH1D((name_+Form("hLamBarF3MCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarF2PCos_[cent][ipt] = new TH1D((name_+Form("hLamBarF2PCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarF3PCos_[cent][ipt] = new TH1D((name_+Form("hLamBarF3PCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaPV2M_[cent][ipt] = new TH1D((name_+Form("hLambdaPV2M_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaMV2P_[cent][ipt] = new TH1D((name_+Form("hLambdaMV2P_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaFV2M_[cent][ipt] = new TH1D((name_+Form("hLambdaFV2M_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaFV2P_[cent][ipt] = new TH1D((name_+Form("hLambdaFV2P_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarPV2M_[cent][ipt] = new TH1D((name_+Form("hLamBarPV2M_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarMV2P_[cent][ipt] = new TH1D((name_+Form("hLamBarMV2P_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarFV2M_[cent][ipt] = new TH1D((name_+Form("hLamBarFV2M_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarFV2P_[cent][ipt] = new TH1D((name_+Form("hLamBarFV2P_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaPV3M_[cent][ipt] = new TH1D((name_+Form("hLambdaPV3M_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaMV3P_[cent][ipt] = new TH1D((name_+Form("hLambdaMV3P_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaFV3M_[cent][ipt] = new TH1D((name_+Form("hLambdaFV3M_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaFV3P_[cent][ipt] = new TH1D((name_+Form("hLambdaFV3P_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarPV3M_[cent][ipt] = new TH1D((name_+Form("hLamBarPV3M_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarMV3P_[cent][ipt] = new TH1D((name_+Form("hLamBarMV3P_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarFV3M_[cent][ipt] = new TH1D((name_+Form("hLamBarFV3M_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarFV3P_[cent][ipt] = new TH1D((name_+Form("hLamBarFV3P_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", int(CentBins[cent]), int(CentBins[cent+1]), pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
            }
        }
    };

    int Fill(const PolarEvent& evt) {
        hCent_->Fill(evt.Cent);
        if ( evt.Cent >= CentBins[NCent] ) return 2;

        int cent = 0;
        while ( int(evt.Cent) >= CentBins[cent+1] ) cent++;

        int HFp2 = bpPb?static_cast<int>(hipPb::HFp2):static_cast<int>(hiPbPb::HFp2);
        int HFm2 = bpPb?static_cast<int>(hipPb::HFm2):static_cast<int>(hiPbPb::HFm2);
        int trackmid2 = bpPb?static_cast<int>(hipPb::trackmid2):static_cast<int>(hiPbPb::trackmid2);

        int HFp3 = bpPb?static_cast<int>(hipPb::HFp3):static_cast<int>(hiPbPb::HFp3);
        int HFm3 = bpPb?static_cast<int>(hipPb::HFm3):static_cast<int>(hiPbPb::HFm3);
        int trackmid3 = bpPb?static_cast<int>(hipPb::trackmid3):static_cast<int>(hiPbPb::trackmid3);

        hHF2pO_[cent]->Fill( (*evt.EPOrg)[HFp2] );
        hHF2mO_[cent]->Fill( (*evt.EPOrg)[HFm2] );
        hHF2pF_[cent]->Fill( (*evt.EPFlat)[HFp2] );
        hHF2mF_[cent]->Fill( (*evt.EPFlat)[HFm2] );

        hHF3pO_[cent]->Fill( (*evt.EPOrg)[HFp3] );
        hHF3mO_[cent]->Fill( (*evt.EPOrg)[HFm3] );
        hHF3pF_[cent]->Fill( (*evt.EPFlat)[HFp3] );
        hHF3mF_[cent]->Fill( (*evt.EPFlat)[HFm3] );

        std::complex rpHFp2   ( (*evt.EPFlatCos)[HFp2],      (*evt.EPFlatSin)[HFp2] );
        std::complex rpHFm2   ( (*evt.EPFlatCos)[HFm2],      (*evt.EPFlatSin)[HFm2] );
        std::complex rpTrkMid2( (*evt.EPFlatCos)[trackmid2], (*evt.EPFlatSin)[trackmid2] );

        std::complex rpHFp3   ( (*evt.EPFlatCos)[HFp3],      (*evt.EPFlatSin)[HFp3] );
        std::complex rpHFm3   ( (*evt.EPFlatCos)[HFm3],      (*evt.EPFlatSin)[HFm3] );
        std::complex rpTrkMid3( (*evt.EPFlatCos)[trackmid3], (*evt.EPFlatSin)[trackmid3] );

        if ( TMath::IsNaN(rpHFp2.real()) or TMath::IsNaN(rpHFp2.imag()) )       {std::cout << " ---> rpHFp2 NaN" << std::endl;    return 1; }
        if ( TMath::IsNaN(rpHFp3.real()) or TMath::IsNaN(rpHFp3.imag()) )       {std::cout << " ---> rpHFp3 NaN" << std::endl;    return 1; }
        if ( TMath::IsNaN(rpHFm2.real()) or TMath::IsNaN(rpHFm2.imag()) )       {std::cout << " ---> rpHFm2 NaN" << std::endl;    return 1; }
        if ( TMath::IsNaN(rpHFm3.real()) or TMath::IsNaN(rpHFm3.imag()) )       {std::cout << " ---> rpHFm3 NaN" << std::endl;    return 1; }
        if ( TMath::IsNaN(rpTrkMid2.real()) or TMath::IsNaN(rpTrkMid2.imag()) ) {std::cout << " ---> rpTrkMid2 NaN" << std::endl; return 1; }
        if ( TMath::IsNaN(rpTrkMid3.real()) or TMath::IsNaN(rpTrkMid3.imag()) ) {std::cout << " ---> rpTrkMid3 NaN" << std::endl; return 1; }

        hReHFp2HFm2_->Fill( evt.Cent, (rpHFp2*std::conj(rpHFm2)).real() );
        hImHFp2HFm2_->Fill( evt.Cent, (rpHFp2*std::conj(rpHFm2)).imag() );
        hReHFp3HFm3_->Fill( evt.Cent, (rpHFp3*std::conj(rpHFm3)).real() );
        hImHFp3HFm3_->Fill( evt.Cent, (rpHFp3*std::conj(rpHFm3)).imag() );

        hReHFp2TrkMid2_->Fill( evt.Cent, (rpHFp2*std::conj(rpTrkMid2)).real() );
        hImHFp2TrkMid2_->Fill( evt.Cent, (rpHFp2*std::conj(rpTrkMid2)).imag() );
        hReHFp3TrkMid3_->Fill( evt.Cent, (rpHFp3*std::conj(rpTrkMid3)).real() );
        hImHFp3TrkMid3_->Fill( evt.Cent, (rpHFp3*std::conj(rpTrkMid3)).imag() );

        hReHFm2TrkMid2_->Fill( evt.Cent, (rpHFm2*std::conj(rpTrkMid2)).real() );
        hImHFm2TrkMid2_->Fill( evt.Cent, (rpHFm2*std::conj(rpTrkMid2)).imag() );
        hReHFm3TrkMid3_->Fill( evt.Cent, (rpHFm3*std::conj(rpTrkMid3)).real() );
        hImHFm3TrkMid3_->Fill( evt.Cent, (rpHFm3*std::conj(rpTrkMid3)).imag() );

        int NLm = 0;
        for ( unsigned int i = 0; i < evt.pdgId->size(); i++ ) {
            if ( ((*evt.pt)[i] < pTbin[0]) or ((*evt.pt)[i] > pTbin[NpT]) ) continue;
            if ( bForward_ ) {
                if ( ( abs((*evt.eta)[i]) < 1.0 ) or ( abs((*evt.eta)[i]) > 2.0 ) ) continue;
            } else {
                if ( ( abs((*evt.eta)[i]) > 1.0 ) ) continue;
            }
            NLm++;
            hLambdaEta_[cent]->Fill( (*evt.eta)[i] );

            int ipt = 0;
            while ( (*evt.pt)[i] > pTbin[ipt+1] ) ipt++;

            double Dphi_m2 = (*evt.phi)[i] - (*evt.EPFlat)[HFm2] ;
            double Dphi_p2 = (*evt.phi)[i] - (*evt.EPFlat)[HFp2] ;
            double Dphi_m3 = (*evt.phi)[i] - (*evt.EPFlat)[HFm3] ;
            double Dphi_p3 = (*evt.phi)[i] - (*evt.EPFlat)[HFp3] ;

            std::complex Q2(TMath::Cos(2.*(*evt.phi)[i]), TMath::Sin(2.*(*evt.phi)[i]));
            std::complex Q3(TMath::Cos(3.*(*evt.phi)[i]), TMath::Sin(3.*(*evt.phi)[i]));

            if ( (*evt.pdgId)[i] == 3122 ) {
                // Lambda
                double cosTheta = (*evt.pPzCM)[i] / sqrt( (*evt.pPzCM)[i]*(*evt.pPzCM)[i] + (*evt.pPxCM)[i]*(*evt.pPxCM)[i] + (*evt.pPyCM)[i]*(*evt.pPyCM)[i] );
                hLambdaMass_[cent][ipt]->Fill( (*evt.mass)[i] );
                if ( (*evt.eta)[i] > 0 ) {
                    // P-side
                    hLambdaP2Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_m2 ) );
                    hLambdaP2Cos2_[cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                    hLambdaP3Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_m3 ) );
                    hLambdaPV2M_  [cent][ipt]->Fill( (*evt.mass)[i], (Q2*std::conj(rpHFm2)).real() );
                    hLambdaPV3M_  [cent][ipt]->Fill( (*evt.mass)[i], (Q3*std::conj(rpHFm3)).real() );
                } else {
                    // M-side
                    hLambdaM2Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_p2 ) );
                    hLambdaM2Cos2_[cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                    hLambdaM3Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_p3 ) );
                    hLambdaMV2P_  [cent][ipt]->Fill( (*evt.mass)[i], (Q2*std::conj(rpHFp2)).real() );
                    hLambdaMV3P_  [cent][ipt]->Fill( (*evt.mass)[i], (Q3*std::conj(rpHFp3)).real() );
                }
                hLambdaF2MCos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_m2 ) );
                hLambdaF2Cos2_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                hLambdaF3MCos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_m3 ) );
                hLambdaF2PCos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_p2 ) );
                hLambdaF3PCos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_p3 ) );
                hLambdaFV2M_   [cent][ipt]->Fill( (*evt.mass)[i], (Q2*std::conj(rpHFm2)).real() );
                hLambdaFV3M_   [cent][ipt]->Fill( (*evt.mass)[i], (Q3*std::conj(rpHFm3)).real() );
                hLambdaFV2P_   [cent][ipt]->Fill( (*evt.mass)[i], (Q2*std::conj(rpHFp2)).real() );
                hLambdaFV3P_   [cent][ipt]->Fill( (*evt.mass)[i], (Q3*std::conj(rpHFp3)).real() );
            }
            if ( (*evt.pdgId)[i] == -3122 ) {
                // anti-Lambda
                double cosTheta = (*evt.nPzCM)[i] / sqrt( (*evt.nPzCM)[i]*(*evt.nPzCM)[i] + (*evt.nPxCM)[i]*(*evt.nPxCM)[i] + (*evt.nPyCM)[i]*(*evt.nPyCM)[i] );
                hLamBarMass_[cent][ipt]->Fill( (*evt.mass)[i] );
                if ( (*evt.eta)[i] > 0 ) {
                    // P-side
                    hLamBarP2Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_m2 ) );
                    hLamBarP2Cos2_[cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                    hLamBarP3Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_m3 ) );
                    hLamBarPV2M_  [cent][ipt]->Fill( (*evt.mass)[i], (Q2*std::conj(rpHFm2)).real() );
                    hLamBarPV3M_  [cent][ipt]->Fill( (*evt.mass)[i], (Q3*std::conj(rpHFm3)).real() );
                } else {
                    // M-side
                    hLamBarM2Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_p2 ) );
                    hLamBarM2Cos2_[cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                    hLamBarM3Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_p3 ) );
                    hLamBarMV2P_  [cent][ipt]->Fill( (*evt.mass)[i], (Q2*std::conj(rpHFp2)).real() );
                    hLamBarMV3P_  [cent][ipt]->Fill( (*evt.mass)[i], (Q3*std::conj(rpHFp3)).real() );
                }
                hLamBarF2MCos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_m2 ) );
                hLamBarF2Cos2_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                hLamBarF3MCos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_m3 ) );
                hLamBarF2PCos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_p2 ) );
                hLamBarF3PCos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_p3 ) );
                hLamBarFV2M_   [cent][ipt]->Fill( (*evt.mass)[i], (Q2*std::conj(rpHFm2)).real() );
                hLamBarFV3M_   [cent][ipt]->Fill( (*evt.mass)[i], (Q3*std::conj(rpHFm3)).real() );
                hLamBarFV2P_   [cent][ipt]->Fill( (*evt.mass)[i], (Q2*std::conj(rpHFp2)).real() );
                hLamBarFV3P_   [cent][ipt]->Fill( (*evt.mass)[i], (Q3*std::conj(rpHFp3)).real() );
            }
        }
        hNLambda_[cent]->Fill(NLm);
        return 0;
    }

    void Write() {
        hCent_->Write("hCent");

        hReHFp2HFm2_   ->Write("hReHFp2HFm2");;
        hReHFp2TrkMid2_->Write("hReHFp2TrkMid2");;
        hReHFm2TrkMid2_->Write("hReHFm2TrkMid2");;
        hImHFp2HFm2_   ->Write("hImHFp2HFm2");;
        hImHFp2TrkMid2_->Write("hImHFp2TrkMid2");;
        hImHFm2TrkMid2_->Write("hImHFm2TrkMid2");;

        hReHFp3HFm3_   ->Write("hReHFp3HFm3");;
        hReHFp3TrkMid3_->Write("hReHFp3TrkMid3");;
        hReHFm3TrkMid3_->Write("hReHFm3TrkMid3");;
        hImHFp3HFm3_   ->Write("hImHFp3HFm3");;
        hImHFp3TrkMid3_->Write("hImHFp3TrkMid3");;
        hImHFm3TrkMid3_->Write("hImHFm3TrkMid3");;

        for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
            hHF2pO_[cent]->Write(Form("hHF2pO_%i", cent));
            hHF2mO_[cent]->Write(Form("hHF2mO_%i", cent));
            hHF2pF_[cent]->Write(Form("hHF2pF_%i", cent));
            hHF2mF_[cent]->Write(Form("hHF2mF_%i", cent));
            hHF3pO_[cent]->Write(Form("hHF3pO_%i", cent));
            hHF3mO_[cent]->Write(Form("hHF3mO_%i", cent));
            hHF3pF_[cent]->Write(Form("hHF3pF_%i", cent));
            hHF3mF_[cent]->Write(Form("hHF3mF_%i", cent));

            hNLambda_[cent]->Write(Form("hNLambda_%i", cent));
            hLambdaEta_[cent]->Write(Form("hLambdaEta_%i", cent));

            for ( int ipt = 0; ipt < NpT; ipt++ ) {
                hLambdaMass_  [cent][ipt]->Write(Form("hLambdaMass_%i_%i", cent, ipt));
                hLamBarMass_  [cent][ipt]->Write(Form("hLamBarMass_%i_%i", cent, ipt));

                hLambdaP2Cos_ [cent][ipt]->Write(Form("hLambdaP2Cos_%i_%i", cent, ipt));
                hLambdaP2Cos2_[cent][ipt]->Write(Form("hLambdaP2Cos2_%i_%i", cent, ipt));
                hLamBarP2Cos_ [cent][ipt]->Write(Form("hLamBarP2Cos_%i_%i ", cent, ipt));
                hLamBarP2Cos2_[cent][ipt]->Write(Form("hLamBarP2Cos2_%i_%i", cent, ipt));
                hLambdaM2Cos_ [cent][ipt]->Write(Form("hLambdaM2Cos_%i_%i", cent, ipt));
                hLambdaM2Cos2_[cent][ipt]->Write(Form("hLambdaM2Cos2_%i_%i", cent, ipt));
                hLamBarM2Cos_ [cent][ipt]->Write(Form("hLamBarM2Cos_%i_%i", cent, ipt));
                hLamBarM2Cos2_[cent][ipt]->Write(Form("hLamBarM2Cos2_%i_%i", cent, ipt));

                hLambdaP3Cos_ [cent][ipt]->Write(Form("hLambdaP3Cos_%i_%i", cent, ipt));
                hLamBarP3Cos_ [cent][ipt]->Write(Form("hLamBarP3Cos_%i_%i", cent, ipt));
                hLambdaM3Cos_ [cent][ipt]->Write(Form("hLambdaM3Cos_%i_%i", cent, ipt));
                hLamBarM3Cos_ [cent][ipt]->Write(Form("hLamBarM3Cos_%i_%i", cent, ipt));

                hLambdaF2MCos_[cent][ipt]->Write(Form("hLambdaF2MCos_%i_%i", cent, ipt));
                hLambdaF2Cos2_[cent][ipt]->Write(Form("hLambdaF2Cos2_%i_%i", cent, ipt));
                hLambdaF3MCos_[cent][ipt]->Write(Form("hLambdaF3MCos_%i_%i", cent, ipt));
                hLambdaF2PCos_[cent][ipt]->Write(Form("hLambdaF2PCos_%i_%i", cent, ipt));
                hLambdaF3PCos_[cent][ipt]->Write(Form("hLambdaF3PCos_%i_%i", cent, ipt));

                hLamBarF2MCos_[cent][ipt]->Write(Form("hLamBarF2MCos_%i_%i", cent, ipt));
                hLamBarF2Cos2_[cent][ipt]->Write(Form("hLamBarF2Cos2_%i_%i", cent, ipt));
                hLamBarF3MCos_[cent][ipt]->Write(Form("hLamBarF3MCos_%i_%i", cent, ipt));
                hLamBarF2PCos_[cent][ipt]->Write(Form("hLamBarF2PCos_%i_%i", cent, ipt));
                hLamBarF3PCos_[cent][ipt]->Write(Form("hLamBarF3PCos_%i_%i", cent, ipt));

                hLambdaPV2M_[cent][ipt]->Write(Form("hLambdaPV2M_%i_%i", cent, ipt));
                hLambdaMV2P_[cent][ipt]->Write(Form("hLambdaMV2P_%i_%i", cent, ipt));
                hLambdaFV2M_[cent][ipt]->Write(Form("hLambdaFV2M_%i_%i", cent, ipt));
                hLambdaFV2P_[cent][ipt]->Write(Form("hLambdaFV2P_%i_%i", cent, ipt));
                hLamBarPV2M_[cent][ipt]->Write(Form("hLamBarPV2M_%i_%i", cent, ipt));
                hLamBarMV2P_[cent][ipt]->Write(Form("hLamBarMV2P_%i_%i", cent, ipt));
                hLamBarFV2M_[cent][ipt]->Write(Form("hLamBarFV2M_%i_%i", cent, ipt));
                hLamBarFV2P_[cent][ipt]->Write(Form("hLamBarFV2P_%i_%i", cent, ipt));
            }
        }
    }

    std::string name_;
    bool    bForward_;

    TH1D*   hCent_;

    vector<TH1D*>   hHF2pO_    ;
    vector<TH1D*>   hHF2mO_    ;
    vector<TH1D*>   hHF2pF_    ;
    vector<TH1D*>   hHF2mF_    ;
    vector<TH1D*>   hHF3pO_    ;
    vector<TH1D*>   hHF3mO_    ;
    vector<TH1D*>   hHF3pF_    ;
    vector<TH1D*>   hHF3mF_    ;
    vector<TH1D*>   hNLambda_  ;
    vector<TH1D*>   hLambdaEta_;

    TH1D*   hReHFp2HFm2_   ;
    TH1D*   hReHFp2TrkMid2_;
    TH1D*   hReHFm2TrkMid2_;
    TH1D*   hImHFp2HFm2_   ;
    TH1D*   hImHFp2TrkMid2_;
    TH1D*   hImHFm2TrkMid2_;

    TH1D*   hReHFp3HFm3_   ;
    TH1D*   hReHFp3TrkMid3_;
    TH1D*   hReHFm3TrkMid3_;
    TH1D*   hImHFp3HFm3_   ;
    TH1D*   hImHFp3TrkMid3_;
    TH1D*   hImHFm3TrkMid3_;

    // Lm mass binning
    vector<vector<TH1D*>>   hLambdaMass_  ;
    vector<vector<TH1D*>>   hLamBarMass_  ;

    vector<vector<TH1D*>>   hLambdaP2Cos_ ;
    vector<vector<TH1D*>>   hLambdaP2Cos2_;
    vector<vector<TH1D*>>   hLamBarP2Cos_ ;
    vector<vector<TH1D*>>   hLamBarP2Cos2_;

    vector<vector<TH1D*>>   hLambdaM2Cos_ ;
    vector<vector<TH1D*>>   hLambdaM2Cos2_;
    vector<vector<TH1D*>>   hLamBarM2Cos_ ;
    vector<vector<TH1D*>>   hLamBarM2Cos2_;

    vector<vector<TH1D*>>   hLambdaP3Cos_ ;
    vector<vector<TH1D*>>   hLamBarP3Cos_ ;
    vector<vector<TH1D*>>   hLambdaM3Cos_ ;
    vector<vector<TH1D*>>   hLamBarM3Cos_ ;

    // full eta range
    vector<vector<TH1D*>>   hLambdaF2MCos_ ;
    vector<vector<TH1D*>>   hLambdaF2Cos2_ ;
    vector<vector<TH1D*>>   hLambdaF3MCos_ ;
    vector<vector<TH1D*>>   hLambdaF2PCos_ ;
    vector<vector<TH1D*>>   hLambdaF3PCos_ ;

    vector<vector<TH1D*>>   hLamBarF2MCos_ ;
    vector<vector<TH1D*>>   hLamBarF2Cos2_ ;
    vector<vector<TH1D*>>   hLamBarF3MCos_ ;
    vector<vector<TH1D*>>   hLamBarF2PCos_ ;
    vector<vector<TH1D*>>   hLamBarF3PCos_ ;

    // v2
    vector<vector<TH1D*>>   hLambdaPV2M_;
    vector<vector<TH1D*>>   hLambdaMV2P_;
    vector<vector<TH1D*>>   hLambdaFV2M_;
    vector<vector<TH1D*>>   hLambdaFV2P_;
    vector<vector<TH1D*>>   hLamBarPV2M_;
    vector<vector<TH1D*>>   hLamBarMV2P_;
    vector<vector<TH1D*>>   hLamBarFV2M_;
    vector<vector<TH1D*>>   hLamBarFV2P_;

    vector<vector<TH1D*>>   hLambdaPV3M_;
    vector<vector<TH1D*>>   hLambdaMV3P_;
    vector<vector<TH1D*>>   hLambdaFV3M_;
    vector<vector<TH1D*>>   hLambdaFV3P_;
    vector<vector<TH1D*>>   hLamBarPV3M_;
    vector<vector<TH1D*>>   hLamBarMV3P_;
    vector<vector<TH1D*>>   hLamBarFV3M_;
    vector<vector<TH1D*>>   hLamBarFV3P_;
};

