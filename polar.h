#include "const.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <vector>
#include <complex>

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

        hReHFp2HFm2_    = new TH1D((name_+"hReHFp2HFm2").c_str(),    "hReHFp2HFm2",    20, 0, 200);
        hReHFp2TrkMid2_ = new TH1D((name_+"hReHFp2TrkMid2").c_str(), "hReHFp2TrkMid2", 20, 0, 200);
        hReHFm2TrkMid2_ = new TH1D((name_+"hReHFm2TrkMid2").c_str(), "hReHFm2TrkMid2", 20, 0, 200);
        hImHFp2HFm2_    = new TH1D((name_+"hImHFp2HFm2").c_str(),    "hImHFp2HFm2",    20, 0, 200);
        hImHFp2TrkMid2_ = new TH1D((name_+"hImHFp2TrkMid2").c_str(), "hImHFp2TrkMid2", 20, 0, 200);
        hImHFm2TrkMid2_ = new TH1D((name_+"hImHFm2TrkMid2").c_str(), "hImHFm2TrkMid2", 20, 0, 200);

        hReHFp3HFm3_    = new TH1D((name_+"hReHFp3HFm3").c_str(),    "hReHFp3HFm3",    20, 0, 200);
        hReHFp3TrkMid3_ = new TH1D((name_+"hReHFp3TrkMid3").c_str(), "hReHFp3TrkMid3", 20, 0, 200);
        hReHFm3TrkMid3_ = new TH1D((name_+"hReHFm3TrkMid3").c_str(), "hReHFm3TrkMid3", 20, 0, 200);
        hImHFp3HFm3_    = new TH1D((name_+"hImHFp3HFm3").c_str(),    "hImHFp3HFm3",    20, 0, 200);
        hImHFp3TrkMid3_ = new TH1D((name_+"hImHFp3TrkMid3").c_str(), "hImHFp3TrkMid3", 20, 0, 200);
        hImHFm3TrkMid3_ = new TH1D((name_+"hImHFm3TrkMid3").c_str(), "hImHFm3TrkMid3", 20, 0, 200);

        for ( int cent = 0; cent < NCentPbPb2018; cent++ ) {
            hHF2pO_[cent] = new TH1D((name_+Form("hHF2pO_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2mO_[cent] = new TH1D((name_+Form("hHF2mO_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2pF_[cent] = new TH1D((name_+Form("hHF2pF_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);
            hHF2mF_[cent] = new TH1D((name_+Form("hHF2mF_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/2., TMath::Pi()/2.);

            hHF3pO_[cent] = new TH1D((name_+Form("hHF3pO_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/3., TMath::Pi()/3.);
            hHF3mO_[cent] = new TH1D((name_+Form("hHF3mO_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/3., TMath::Pi()/3.);
            hHF3pF_[cent] = new TH1D((name_+Form("hHF3pF_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/3., TMath::Pi()/3.);
            hHF3mF_[cent] = new TH1D((name_+Form("hHF3mF_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 72, -TMath::Pi()/3., TMath::Pi()/3.);

            hNLambda_[cent] = new TH1D((name_+Form("hNLambda_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 20, 0, 20);
            hLambdaEta_[cent] = new TH1D((name_+Form("hLambdaEta_%i", cent)).c_str(), Form("CentBin %i-%i", CentPbPb2018[cent], CentPbPb2018[cent+1]), 50, -2.5, 2.5);

            for ( int ipt = 0; ipt < NpT; ipt++ ) {
                hLambdaMass_[cent][ipt] = new TH1D((name_+Form("hLambdaMass_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarMass_[cent][ipt] = new TH1D((name_+Form("hLamBarMass_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaP2Cos_ [cent][ipt] = new TH1D((name_+Form("hLambdaP2Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaP2Cos2_[cent][ipt] = new TH1D((name_+Form("hLambdaP2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarP2Cos_ [cent][ipt] = new TH1D((name_+Form("hLamBarP2Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarP2Cos2_[cent][ipt] = new TH1D((name_+Form("hLamBarP2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaM2Cos_ [cent][ipt] = new TH1D((name_+Form("hLambdaM2Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLambdaM2Cos2_[cent][ipt] = new TH1D((name_+Form("hLambdaM2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarM2Cos_ [cent][ipt] = new TH1D((name_+Form("hLamBarM2Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarM2Cos2_[cent][ipt] = new TH1D((name_+Form("hLamBarM2Cos2_%i_%i", cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaP3Cos_ [cent][ipt] = new TH1D((name_+Form("hLambdaP3Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarP3Cos_ [cent][ipt] = new TH1D((name_+Form("hLamBarP3Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);

                hLambdaM3Cos_ [cent][ipt] = new TH1D((name_+Form("hLambdaM3Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
                hLamBarM3Cos_ [cent][ipt] = new TH1D((name_+Form("hLamBarM3Cos_%i_%i",  cent, ipt)).c_str(), Form("CentBin %i-%i, pT %.1f-%.1f", CentPbPb2018[cent], CentPbPb2018[cent+1], pTbin[ipt], pTbin[ipt+1]), N_LmBins, Lm_mass_min, Lm_mass_max);
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

        hHF3pO_[cent]->Fill( (*evt.EPOrg)[hi::HFp3] );
        hHF3mO_[cent]->Fill( (*evt.EPOrg)[hi::HFm3] );
        hHF3pF_[cent]->Fill( (*evt.EPFlat)[hi::HFp3] );
        hHF3mF_[cent]->Fill( (*evt.EPFlat)[hi::HFm3] );

        std::complex rpHFp2   ( (*evt.EPFlatCos)[hi::HFp2],      (*evt.EPFlatSin)[hi::HFp2] );
        std::complex rpHFm2   ( (*evt.EPFlatCos)[hi::HFm2],      (*evt.EPFlatSin)[hi::HFm2] );
        std::complex rpTrkMid2( (*evt.EPFlatCos)[hi::trackmid2], (*evt.EPFlatSin)[hi::trackmid2] );

        std::complex rpHFp3   ( (*evt.EPFlatCos)[hi::HFp3],      (*evt.EPFlatSin)[hi::HFp3] );
        std::complex rpHFm3   ( (*evt.EPFlatCos)[hi::HFm3],      (*evt.EPFlatSin)[hi::HFm3] );
        std::complex rpTrkMid3( (*evt.EPFlatCos)[hi::trackmid3], (*evt.EPFlatSin)[hi::trackmid3] );

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

            double Dphi_m2 = (*evt.phi)[i] - (*evt.EPFlat)[hi::HFm2] ;
            double Dphi_p2 = (*evt.phi)[i] - (*evt.EPFlat)[hi::HFp2] ;
            double Dphi_m3 = (*evt.phi)[i] - (*evt.EPFlat)[hi::HFm3] ;
            double Dphi_p3 = (*evt.phi)[i] - (*evt.EPFlat)[hi::HFp3] ;

            if ( (*evt.pdgId)[i] == 3122 ) {
                // Lambda
                double cosTheta = (*evt.pPzCM)[i] / sqrt( (*evt.pPzCM)[i]*(*evt.pPzCM)[i] + (*evt.pPxCM)[i]*(*evt.pPxCM)[i] + (*evt.pPyCM)[i]*(*evt.pPyCM)[i] );
                hLambdaMass_[cent][ipt]->Fill( (*evt.mass)[i] );
                if ( (*evt.eta)[i] > 0 ) {
                    // P-side
                    hLambdaP2Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_m2 ) );
                    hLambdaP2Cos2_[cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                    hLambdaP3Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_m3 ) );
                } else {
                    // M-side
                    hLambdaM2Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_p2 ) );
                    hLambdaM2Cos2_[cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                    hLambdaM3Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_p3 ) );
                }
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
                } else {
                    // M-side
                    hLamBarM2Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 2*Dphi_p2 ) );
                    hLamBarM2Cos2_[cent][ipt]->Fill( (*evt.mass)[i], cosTheta*cosTheta );
                    hLamBarM3Cos_ [cent][ipt]->Fill( (*evt.mass)[i], cosTheta*TMath::Sin( 3*Dphi_p3 ) );
                }
            }
        }
        hNLambda_[cent]->Fill(NLm);
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

        for ( int cent = 0; cent < NCentPbPb2018; cent++ ) {
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
    TH1D*   hHF3pO_[NCentPbPb2018];
    TH1D*   hHF3mO_[NCentPbPb2018];
    TH1D*   hHF3pF_[NCentPbPb2018];
    TH1D*   hHF3mF_[NCentPbPb2018];

    TH1D*   hNLambda_[NCentPbPb2018];
    TH1D*   hLambdaEta_[NCentPbPb2018];

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
    TH1D*   hLambdaMass_[NCentPbPb2018][NpT];
    TH1D*   hLamBarMass_[NCentPbPb2018][NpT];

    TH1D*   hLambdaP2Cos_ [NCentPbPb2018][NpT];
    TH1D*   hLambdaP2Cos2_[NCentPbPb2018][NpT];
    TH1D*   hLamBarP2Cos_ [NCentPbPb2018][NpT];
    TH1D*   hLamBarP2Cos2_[NCentPbPb2018][NpT];
    TH1D*   hLambdaM2Cos_ [NCentPbPb2018][NpT];
    TH1D*   hLambdaM2Cos2_[NCentPbPb2018][NpT];
    TH1D*   hLamBarM2Cos_ [NCentPbPb2018][NpT];
    TH1D*   hLamBarM2Cos2_[NCentPbPb2018][NpT];

    TH1D*   hLambdaP3Cos_ [NCentPbPb2018][NpT];
    TH1D*   hLamBarP3Cos_ [NCentPbPb2018][NpT];
    TH1D*   hLambdaM3Cos_ [NCentPbPb2018][NpT];
    TH1D*   hLamBarM3Cos_ [NCentPbPb2018][NpT];
};

