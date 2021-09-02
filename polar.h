#include "const.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <vector>

using namespace std;
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
    PolarData(std::string name, bool bFwd = false) {
        name_ = name;
        bForward_ = bFwd;
        hCent_ = new TH1D((name_+"hCent").c_str(), "hCent", 200, 0, 200);
        for ( int cent = 0; cent < NCentPbPb2018; cent++ ) {
            hHF2pO_[cent] = new TH1D((name_+Form("hHF2pO_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2mO_[cent] = new TH1D((name_+Form("hHF2mO_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2pF_[cent] = new TH1D((name_+Form("hHF2pF_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2mF_[cent] = new TH1D((name_+Form("hHF2mF_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2Res_[cent] = new TH1D((name_+Form("hHF2Res_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);

            hNLambda_[cent] = new TH1D((name_+Form("hNLambda_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 20, 0, 20);

            for ( int ipt = 0; ipt < NpT; ipt++ ) {
                h2LambdaPCos_[cent][ipt] = new TH2D((name_+Form("h2LambdaPCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], ptBin[ipt], ptBin[ipt+1]), 12, -TMath::Pi(), TMath::Pi(), 20, -1., 1. );
                h2LambdaMCos_[cent][ipt] = new TH2D((name_+Form("h2LambdaMCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], ptBin[ipt], ptBin[ipt+1]), 12, -TMath::Pi(), TMath::Pi(), 20, -1., 1. );
                h2LamBarPCos_[cent][ipt] = new TH2D((name_+Form("h2LamBarPCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], ptBin[ipt], ptBin[ipt+1]), 12, -TMath::Pi(), TMath::Pi(), 20, -1., 1. );
                h2LamBarMCos_[cent][ipt] = new TH2D((name_+Form("h2LamBarMCos_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], ptBin[ipt], ptBin[ipt+1]), 12, -TMath::Pi(), TMath::Pi(), 20, -1., 1. );

                hLambdaPhiP_[cent][ipt] = new TH1D((name_+Form("hLambdaPhiP_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], ptBin[ipt], ptBin[ipt+1]), 72, -TMath::Pi(), TMath::Pi());
                hLamBarPhiP_[cent][ipt] = new TH1D((name_+Form("hLamBarPhiP_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], ptBin[ipt], ptBin[ipt+1]), 72, -TMath::Pi(), TMath::Pi());
                hLambdaPhiM_[cent][ipt] = new TH1D((name_+Form("hLambdaPhiM_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], ptBin[ipt], ptBin[ipt+1]), 72, -TMath::Pi(), TMath::Pi());
                hLamBarPhiM_[cent][ipt] = new TH1D((name_+Form("hLamBarPhiM_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], ptBin[ipt], ptBin[ipt+1]), 72, -TMath::Pi(), TMath::Pi());
            }
        }
    };

    void Fill(const PolarEvent& evt) {
        hCent_->Fill(evt.Cent);
        if ( evt.Cent >= CentPbPb2018[NCentPbPb2018] ) return;

        int cent = 0;
        while ( int(evt.Cent) >= CentPbPb2018[cent+1] ) cent++;

        hHF2pO_[cent]->Fill( (*evt.EPOrg)[hi::HFp2] );
        hHF2mO_[cent]->Fill( (*evt.EPOrg)[hi::HFm2] );
        hHF2pF_[cent]->Fill( (*evt.EPFlat)[hi::HFp2] );
        hHF2mF_[cent]->Fill( (*evt.EPFlat)[hi::HFm2] );
        hHF2Res_[cent]->Fill( (*evt.EPFlat)[hi::HFp2] - (*evt.EPFlat)[hi::HFm2] );

        int NLm = 0;
        for ( unsigned int i = 0; i < evt.pdgId->size(); i++ ) {
            if ( ((*evt.pt)[i] < ptBin[0]) or ((*evt.pt)[i] > ptBin[NpT]) ) continue;
            if ( bForward_ ) {
                if ( ( abs((*evt.eta)[i]) < 1.0 ) or ( abs((*evt.eta)[i]) > 2.0 ) ) continue;
            } else {
                if ( ( abs((*evt.eta)[i]) > 1.0 ) ) continue;
            }
            NLm++;
            int ipt = 0;
            while ( (*evt.pt)[i] > ptBin[ipt+1] ) ipt++;

            if ( (*evt.eta)[i] > 0 ) {
                // pos eta
                if ( (*evt.pdgId)[i] == 3122 ) {
                    // Lambda
                    double cosTheta = (*evt.pPzCM)[i] / sqrt( (*evt.pPzCM)[i]*(*evt.pPzCM)[i] + (*evt.pPxCM)[i]*(*evt.pPxCM)[i] + (*evt.pPyCM)[i]*(*evt.pPyCM)[i] );
                    double Dphi = roundPi( (*evt.phi)[i] - (*evt.EPFlat)[hi::HFm2] );
                    h2LambdaPCos_[cent][ipt]->Fill( Dphi, cosTheta );
                    hLambdaPhiP_[cent][ipt]->Fill(Dphi);
                }
                if ( (*evt.pdgId)[i] == -3122 ) {
                    // anti-Lambda
                    double cosTheta = (*evt.nPzCM)[i] / sqrt( (*evt.nPzCM)[i]*(*evt.nPzCM)[i] + (*evt.nPxCM)[i]*(*evt.nPxCM)[i] + (*evt.nPyCM)[i]*(*evt.nPyCM)[i] );
                    double Dphi = roundPi( (*evt.phi)[i] - (*evt.EPFlat)[hi::HFm2] );
                    h2LamBarPCos_[cent][ipt]->Fill( Dphi, cosTheta );
                    hLamBarPhiP_[cent][ipt]->Fill(Dphi);
                }
            } else {
                // neg eta
                if ( (*evt.pdgId)[i] == 3122 ) {
                    // Lambda
                    double cosTheta = (*evt.pPzCM)[i] / sqrt( (*evt.pPzCM)[i]*(*evt.pPzCM)[i] + (*evt.pPxCM)[i]*(*evt.pPxCM)[i] + (*evt.pPyCM)[i]*(*evt.pPyCM)[i] );
                    double Dphi = roundPi( (*evt.phi)[i] - (*evt.EPFlat)[hi::HFp2] );
                    h2LambdaMCos_[cent][ipt]->Fill( Dphi, cosTheta );
                    hLambdaPhiM_[cent][ipt]->Fill(Dphi);
                }
                if ( (*evt.pdgId)[i] == -3122 ) {
                    // anti-Lambda
                    double cosTheta = (*evt.nPzCM)[i] / sqrt( (*evt.nPzCM)[i]*(*evt.nPzCM)[i] + (*evt.nPxCM)[i]*(*evt.nPxCM)[i] + (*evt.nPyCM)[i]*(*evt.nPyCM)[i] );
                    double Dphi = roundPi( (*evt.phi)[i] - (*evt.EPFlat)[hi::HFp2] );
                    h2LamBarMCos_[cent][ipt]->Fill( Dphi, cosTheta );
                    hLamBarPhiM_[cent][ipt]->Fill(Dphi);
                }
            }
        }
        hNLambda_[cent]->Fill(NLm);
    }

    void Write() {
        hCent_->Write("hCent");

        for ( int cent = 0; cent < NCentPbPb2018; cent++ ) {
            hHF2pO_[cent]->Write(Form("hHF2pO_%i", cent));
            hHF2mO_[cent]->Write(Form("hHF2mO_%i", cent));
            hHF2pF_[cent]->Write(Form("hHF2pF_%i", cent));
            hHF2mF_[cent]->Write(Form("hHF2mF_%i", cent));
            hHF2Res_[cent]->Write(Form("hHF2Res_%i", cent));

            hNLambda_[cent]->Write(Form("hNLambda_%i", cent));

            for ( int ipt = 0; ipt < NpT; ipt++ ) {
                h2LambdaPCos_[cent][ipt]->Write(Form("h2LambdaPCos_%i_%i", cent, ipt));
                h2LamBarPCos_[cent][ipt]->Write(Form("h2LambdaMCos_%i_%i", cent, ipt));
                h2LambdaMCos_[cent][ipt]->Write(Form("h2LamBarPCos_%i_%i", cent, ipt));
                h2LamBarMCos_[cent][ipt]->Write(Form("h2LamBarMCos_%i_%i", cent, ipt));

                hLambdaPhiP_[cent][ipt]->Write(Form("hLambdaPhiP_%i_%i", cent, ipt));
                hLamBarPhiP_[cent][ipt]->Write(Form("hLamBarPhiP_%i_%i", cent, ipt));
                hLambdaPhiM_[cent][ipt]->Write(Form("hLambdaPhiM_%i_%i", cent, ipt));
                hLamBarPhiM_[cent][ipt]->Write(Form("hLamBarPhiM_%i_%i", cent, ipt));
            }
        }
    }

    std::string name_;
    bool    bForward_;
    TH1D*   hCent_;

    TH1D*   hHF2pO_[NCentPbPb2018];
    TH1D*   hHF2mO_[NCentPbPb2018];
    TH1D*   hHF2pF_[NCentPbPb2018];
    TH1D*   hHF2mF_[NCentPbPb2018];
    TH1D*   hHF2Res_[NCentPbPb2018];

    TH1D*   hNLambda_[NCentPbPb2018];

    TH2D*   h2LambdaPCos_[NCentPbPb2018][NpT];
    TH2D*   h2LamBarPCos_[NCentPbPb2018][NpT];
    TH2D*   h2LambdaMCos_[NCentPbPb2018][NpT];
    TH2D*   h2LamBarMCos_[NCentPbPb2018][NpT];

    TH1D*   hLambdaPhiP_[NCentPbPb2018][NpT];
    TH1D*   hLamBarPhiP_[NCentPbPb2018][NpT];
    TH1D*   hLambdaPhiM_[NCentPbPb2018][NpT];
    TH1D*   hLamBarPhiM_[NCentPbPb2018][NpT];
};

