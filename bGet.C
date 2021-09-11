#include "label.h"
#include "const.h"

#include <TH1D.h>
#include <TFile.h>

void bGet(int s1 = 0, int s3 = 10)
{
    std::cout << " s1 = " << s1 << " s3 = " << s3 << std::endl;
    TH1::SetDefaultSumw2();

    TFile* f = new TFile(Form("%s/output_%i.root", ftxt[s1], s3));

    //
    TH1D*   hCent_ = 0;
    TH1D*   hCent2_ = 0;

    TH1D*   hHF2pO_[NCentPbPb2018] = {};
    TH1D*   hHF2mO_[NCentPbPb2018] = {};
    TH1D*   hHF2pF_[NCentPbPb2018] = {};
    TH1D*   hHF2mF_[NCentPbPb2018] = {};
    TH1D*   hHF3pO_[NCentPbPb2018] = {};
    TH1D*   hHF3mO_[NCentPbPb2018] = {};
    TH1D*   hHF3pF_[NCentPbPb2018] = {};
    TH1D*   hHF3mF_[NCentPbPb2018] = {};

    TH1D*   hNLambda_[NCentPbPb2018] = {};
    TH1D*   hLambdaEta_[NCentPbPb2018] = {};

    TH1D*   hReHFp2HFm2_    = 0;
    TH1D*   hReHFp2TrkMid2_ = 0;
    TH1D*   hReHFm2TrkMid2_ = 0;
    TH1D*   hImHFp2HFm2_    = 0;
    TH1D*   hImHFp2TrkMid2_ = 0;
    TH1D*   hImHFm2TrkMid2_ = 0;

    TH1D*   hReHFp3HFm3_    = 0;
    TH1D*   hReHFp3TrkMid3_ = 0;
    TH1D*   hReHFm3TrkMid3_ = 0;
    TH1D*   hImHFp3HFm3_    = 0;
    TH1D*   hImHFp3TrkMid3_ = 0;
    TH1D*   hImHFm3TrkMid3_ = 0;

    // Lm mass binning
    TH1D*   hLambdaMass_[NCentPbPb2018][NpT] = {};
    TH1D*   hLamBarMass_[NCentPbPb2018][NpT] = {};

    TH1D*   hLambdaP2Cos_ [NCentPbPb2018][NpT] = {};
    TH1D*   hLambdaP2Cos2_[NCentPbPb2018][NpT] = {};
    TH1D*   hLamBarP2Cos_ [NCentPbPb2018][NpT] = {};
    TH1D*   hLamBarP2Cos2_[NCentPbPb2018][NpT] = {};
    TH1D*   hLambdaM2Cos_ [NCentPbPb2018][NpT] = {};
    TH1D*   hLambdaM2Cos2_[NCentPbPb2018][NpT] = {};
    TH1D*   hLamBarM2Cos_ [NCentPbPb2018][NpT] = {};
    TH1D*   hLamBarM2Cos2_[NCentPbPb2018][NpT] = {};

    TH1D*   hLambdaP3Cos_ [NCentPbPb2018][NpT] = {};
    TH1D*   hLamBarP3Cos_ [NCentPbPb2018][NpT] = {};
    TH1D*   hLambdaM3Cos_ [NCentPbPb2018][NpT] = {};
    TH1D*   hLamBarM3Cos_ [NCentPbPb2018][NpT] = {};

    ////
    hCent_ = (TH1D*) f->Get(Form("output_%i_%i/hCent", s3, s3));
    hCent2_ = (TH1D*) hCent_->Rebin(10, "hCent2");

    hReHFp2HFm2_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2HFm2_", s3, s3));
    hReHFp2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2TrkMid2_", s3, s3));
    hReHFm2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm2TrkMid2_", s3, s3));
    hImHFp2HFm2_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2HFm2_", s3, s3));
    hImHFp2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2TrkMid2_", s3, s3));
    hImHFm2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm2TrkMid2_", s3, s3));

    hReHFp3HFm3_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3HFm3_", s3, s3));
    hReHFp3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3TrkMid3_", s3, s3));
    hReHFm3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm3TrkMid3_", s3, s3));
    hImHFp3HFm3_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3HFm3_", s3, s3));
    hImHFp3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3TrkMid3_", s3, s3));
    hImHFm3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm3TrkMid3_", s3, s3));

    for ( int cent = 0; cent < NCentPbPb2018; cent++ ) {

    hHF2pO_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hHF2pO_%i", s3, s3, cent));
    hHF2mO_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hHF2mO_%i", s3, s3, cent));
    hHF2pF_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hHF2pF_%i", s3, s3, cent));
    hHF2mF_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hHF2mF_%i", s3, s3, cent));
    hHF3pO_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hHF3pO_%i", s3, s3, cent));
    hHF3mO_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hHF3mO_%i", s3, s3, cent));
    hHF3pF_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hHF3pF_%i", s3, s3, cent));
    hHF3mF_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hHF3mF_%i", s3, s3, cent));

    hNLambda_  [cent] = (TH1D*) f->Get(Form("output_%i_%i/hNLambda_%i", s3, s3, cent));
    hLambdaEta_[cent] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaEta_%i", s3, s3, cent));

        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaMass_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaMass_%i_%i", s3, s3, cent, ipt));
            hLamBarMass_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarMass_%i_%i", s3, s3, cent, ipt));

            hLambdaP2Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaP2Cos_%i_%i", s3, s3, cent, ipt));
            hLambdaP2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaP2Cos2_%i_%i", s3, s3, cent, ipt));
            hLamBarP2Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarP2Cos_%i_%i", s3, s3, cent, ipt));
            hLamBarP2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarP2Cos2_%i_%i", s3, s3, cent, ipt));
            hLambdaM2Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaM2Cos_%i_%i", s3, s3, cent, ipt));
            hLambdaM2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaM2Cos2_%i_%i", s3, s3, cent, ipt));
            hLamBarM2Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarM2Cos_%i_%i", s3, s3, cent, ipt));
            hLamBarM2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarM2Cos2_%i_%i", s3, s3, cent, ipt));

            hLambdaP3Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaP3Cos_%i_%i", s3, s3, cent, ipt));
            hLamBarP3Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarP3Cos_%i_%i", s3, s3, cent, ipt));
            hLambdaM3Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaM3Cos_%i_%i", s3, s3, cent, ipt));
            hLamBarM3Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarM3Cos_%i_%i", s3, s3, cent, ipt));
        }
    }
    ////
    //
    // EP
    hReHFp2HFm2_   ->Divide(hCent2_);
    hReHFp2TrkMid2_->Divide(hCent2_);
    hReHFm2TrkMid2_->Divide(hCent2_);
    hImHFp2HFm2_   ->Divide(hCent2_);
    hImHFp2TrkMid2_->Divide(hCent2_);
    hImHFm2TrkMid2_->Divide(hCent2_);

    hReHFp3HFm3_   ->Divide(hCent2_);
    hReHFp3TrkMid3_->Divide(hCent2_);
    hReHFm3TrkMid3_->Divide(hCent2_);
    hImHFp3HFm3_   ->Divide(hCent2_);
    hImHFp3TrkMid3_->Divide(hCent2_);
    hImHFm3TrkMid3_->Divide(hCent2_);

    // Cos2
    for ( int cent = 0; cent < NCentPbPb2018; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaP2Cos2_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBarP2Cos2_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLambdaM2Cos2_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBarM2Cos2_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
        }
    }


    // Save
    TFile * fsave = new TFile(Form("%s/bGet_%i.root", ftxt[s1], s3), "recreate");
    hCent_         ->Write("hCent");
    hCent2_        ->Write("hCent2");
    hReHFp2HFm2_   ->Write("hReHFp2HFm2");
    hReHFp2TrkMid2_->Write("hReHFp2TrkMid2");
    hReHFm2TrkMid2_->Write("hReHFm2TrkMid2");
    hImHFp2HFm2_   ->Write("hImHFp2HFm2");
    hImHFp2TrkMid2_->Write("hImHFp2TrkMid2");
    hImHFm2TrkMid2_->Write("hImHFm2TrkMid2");

    hReHFp3HFm3_   ->Write("hReHFp3HFm3");
    hReHFp3TrkMid3_->Write("hReHFp3TrkMid3");
    hReHFm3TrkMid3_->Write("hReHFm3TrkMid3");
    hImHFp3HFm3_   ->Write("hImHFp3HFm3");
    hImHFp3TrkMid3_->Write("hImHFp3TrkMid3");
    hImHFm3TrkMid3_->Write("hImHFm3TrkMid3");

    for ( int cent = 0; cent < NCentPbPb2018; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaP2Cos2_[cent][ipt]->Write(Form("hLambdaP2Cos2_%i_%i", cent, ipt));
            hLamBarP2Cos2_[cent][ipt]->Write(Form("hLamBarP2Cos2_%i_%i", cent, ipt));
            hLambdaM2Cos2_[cent][ipt]->Write(Form("hLambdaM2Cos2_%i_%i", cent, ipt));
            hLamBarM2Cos2_[cent][ipt]->Write(Form("hLamBarM2Cos2_%i_%i", cent, ipt));
        }
    }
}
