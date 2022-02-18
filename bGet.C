#include "label.h"
#include "const.h"
#include "complex"

#include <TH1D.h>
#include <TFile.h>

bool bDebug = false;

void bGet(int s1 = 0, int s3 = 10)
{
    std::cout << " s1 = " << s1 << " s3 = " << s3 << std::endl;
    TH1::SetDefaultSumw2();
    setFlags(s1);

    int NCent = 0;
    const double* CentBins = nullptr;
    if ( bPbPb ) {
        NCent = NCentPbPb2018;
        CentBins = CentPbPb2018;
    } else if ( bpPb ) {
        NCent = NCentpPb2016;
        CentBins = CentNoffCutPA8TeV4;
    }
    if ( bPbPb ) {
        cout << "  --> PbPb" << endl;
    }
    if ( bpPb ) {
        cout << "  --> pPb";
        if ( bpPbReverse ) {
            cout << " reverse" << endl;
        } else {
            cout << endl;
        }
    }

    TFile* f = new TFile(Form("%s/output_%i.root", ftxt[s1], s3));

    //
    TH1D*   hCent_ = nullptr;
    TH1D*   hCent2_ = nullptr;

    // EP
    // O - Original
    // F - Flatten
    vector<TH1D*>   hHF2pO_    (NCent);
    vector<TH1D*>   hHF2mO_    (NCent);
    vector<TH1D*>   hHF2pF_    (NCent);
    vector<TH1D*>   hHF2mF_    (NCent);
    vector<TH1D*>   hHF3pO_    (NCent);
    vector<TH1D*>   hHF3mO_    (NCent);
    vector<TH1D*>   hHF3pF_    (NCent);
    vector<TH1D*>   hHF3mF_    (NCent);

    vector<TH1D*>   hNLambda_  (NCent);
    vector<TH1D*>   hLambdaEta_(NCent);

    // EP Resolution Parts
    // Re - real
    // Im - imag
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

    TH1D*   hHF2pRes_ = 0;
    TH1D*   hHF2mRes_ = 0;
    TH1D*   hHF3pRes_ = 0;
    TH1D*   hHF3mRes_ = 0;

    // SP Resolution Parts
    // Re - real
    // Im - imag
    TH1D*   hReHFp2HFm2SP_    = 0;
    TH1D*   hReHFp2TrkMid2SP_ = 0;
    TH1D*   hReHFm2TrkMid2SP_ = 0;
    TH1D*   hImHFp2HFm2SP_    = 0;
    TH1D*   hImHFp2TrkMid2SP_ = 0;
    TH1D*   hImHFm2TrkMid2SP_ = 0;

    TH1D*   hReHFp3HFm3SP_    = 0;
    TH1D*   hReHFp3TrkMid3SP_ = 0;
    TH1D*   hReHFm3TrkMid3SP_ = 0;
    TH1D*   hImHFp3HFm3SP_    = 0;
    TH1D*   hImHFp3TrkMid3SP_ = 0;
    TH1D*   hImHFm3TrkMid3SP_ = 0;

    TH1D*   hHF2pResSP_ = 0;
    TH1D*   hHF2mResSP_ = 0;
    TH1D*   hHF3pResSP_ = 0;
    TH1D*   hHF3mResSP_ = 0;

    // Lm mass binning
    // [cent][ipt]
    vector<vector<TH1D*>>   hLambdaMass_   (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMass_   (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMassP_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMassP_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMassM_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMassM_  (NCent, vector<TH1D*>(NpT));

    // Polar P(M)-side wrt HFm(HFp)
    vector<vector<TH1D*>>   hLambdaP2Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2Cos_  (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaP3Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM3Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP3Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM3Cos_  (NCent, vector<TH1D*>(NpT));

    // full eta
    vector<vector<TH1D*>>   hLambdaF2MCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3MCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2PCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3PCos_ (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLamBarF2MCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3MCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2PCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3PCos_ (NCent, vector<TH1D*>(NpT));

    // Polar cosTheta^2 -- acceptance term
    vector<vector<TH1D*>>   hLambdaP2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2Cos2_ (NCent, vector<TH1D*>(NpT));

    // Polar SP
    vector<vector<TH1D*>>   hLambdaP2CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaP3CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM3CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2MCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3MCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2PCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3PCosSP_(NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLamBarP2CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP3CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM3CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2MCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3MCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2PCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3PCosSP_(NCent, vector<TH1D*>(NpT));

    // v2
    vector<vector<TH1D*>>   hLambdaPV2M_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMV2P_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV2M_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV2P_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarPV2M_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMV2P_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV2M_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV2P_ (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaPV3M_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMV3P_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV3M_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV3P_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarPV3M_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMV3P_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV3M_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV3P_ (NCent, vector<TH1D*>(NpT));
    ////
    hCent_ = (TH1D*) f->Get(Form("output_%i_%i/hCent", s3, s3));
    hCent2_ = (TH1D*) hCent_->Rebin(NCent, "hCent2", CentBins);

    hReHFp2HFm2_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2HFm2", s3, s3));
    hReHFp2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2TrkMid2", s3, s3));
    hReHFm2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm2TrkMid2", s3, s3));
    hImHFp2HFm2_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2HFm2", s3, s3));
    hImHFp2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2TrkMid2", s3, s3));
    hImHFm2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm2TrkMid2", s3, s3));

    hReHFp3HFm3_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3HFm3", s3, s3));
    hReHFp3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3TrkMid3", s3, s3));
    hReHFm3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm3TrkMid3", s3, s3));
    hImHFp3HFm3_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3HFm3", s3, s3));
    hImHFp3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3TrkMid3", s3, s3));
    hImHFm3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm3TrkMid3", s3, s3));

    hReHFp2HFm2SP_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2HFm2SP", s3, s3));
    hReHFp2TrkMid2SP_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2TrkMid2SP", s3, s3));
    hReHFm2TrkMid2SP_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm2TrkMid2SP", s3, s3));
    hImHFp2HFm2SP_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2HFm2SP", s3, s3));
    hImHFp2TrkMid2SP_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2TrkMid2SP", s3, s3));
    hImHFm2TrkMid2SP_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm2TrkMid2SP", s3, s3));

    hReHFp3HFm3SP_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3HFm3SP", s3, s3));
    hReHFp3TrkMid3SP_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3TrkMid3SP", s3, s3));
    hReHFm3TrkMid3SP_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm3TrkMid3SP", s3, s3));
    hImHFp3HFm3SP_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3HFm3SP", s3, s3));
    hImHFp3TrkMid3SP_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3TrkMid3SP", s3, s3));
    hImHFm3TrkMid3SP_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm3TrkMid3SP", s3, s3));

    hHF2pRes_ = new TH1D("hHF2pRes_", "hHF2pRes_", NCent, CentBins);
    hHF2mRes_ = new TH1D("hHF2mRes_", "hHF2mRes_", NCent, CentBins);
    hHF3pRes_ = new TH1D("hHF3pRes_", "hHF3pRes_", NCent, CentBins);
    hHF3mRes_ = new TH1D("hHF3mRes_", "hHF3mRes_", NCent, CentBins);

    hHF2pResSP_ = new TH1D("hHF2pResSP_", "hHF2pResSP_", NCent, CentBins);
    hHF2mResSP_ = new TH1D("hHF2mResSP_", "hHF2mResSP_", NCent, CentBins);
    hHF3pResSP_ = new TH1D("hHF3pResSP_", "hHF3pResSP_", NCent, CentBins);
    hHF3mResSP_ = new TH1D("hHF3mResSP_", "hHF3mResSP_", NCent, CentBins);

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
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
            hLambdaMassP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaMassP_%i_%i", s3, s3, cent, ipt));
            hLamBarMassP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarMassP_%i_%i", s3, s3, cent, ipt));
            hLambdaMassM_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaMassM_%i_%i", s3, s3, cent, ipt));
            hLamBarMassM_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarMassM_%i_%i", s3, s3, cent, ipt));

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

            hLambdaP2CosSP_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaP2CosSP_%i_%i", s3, s3, cent, ipt));
            hLambdaP3CosSP_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaP3CosSP_%i_%i", s3, s3, cent, ipt));
            hLambdaM2CosSP_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaM2CosSP_%i_%i", s3, s3, cent, ipt));
            hLambdaM3CosSP_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaM3CosSP_%i_%i", s3, s3, cent, ipt));
            hLambdaF2MCosSP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF2MCosSP_%i_%i", s3, s3, cent, ipt));
            hLambdaF3MCosSP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF3MCosSP_%i_%i", s3, s3, cent, ipt));
            hLambdaF2PCosSP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF2PCosSP_%i_%i", s3, s3, cent, ipt));
            hLambdaF3PCosSP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF3PCosSP_%i_%i", s3, s3, cent, ipt));

            hLamBarP2CosSP_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarP2CosSP_%i_%i", s3, s3, cent, ipt));
            hLamBarP3CosSP_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarP3CosSP_%i_%i", s3, s3, cent, ipt));
            hLamBarM2CosSP_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarM2CosSP_%i_%i", s3, s3, cent, ipt));
            hLamBarM3CosSP_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarM3CosSP_%i_%i", s3, s3, cent, ipt));
            hLamBarF2MCosSP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF2MCosSP_%i_%i", s3, s3, cent, ipt));
            hLamBarF3MCosSP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF3MCosSP_%i_%i", s3, s3, cent, ipt));
            hLamBarF2PCosSP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF2PCosSP_%i_%i", s3, s3, cent, ipt));
            hLamBarF3PCosSP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF3PCosSP_%i_%i", s3, s3, cent, ipt));

            hLambdaF2MCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF2MCos_%i_%i", s3, s3, cent, ipt));
            hLambdaF2Cos2_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF2Cos2_%i_%i", s3, s3, cent, ipt));
            hLambdaF3MCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF3MCos_%i_%i", s3, s3, cent, ipt));
            hLambdaF2PCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF2PCos_%i_%i", s3, s3, cent, ipt));
            hLambdaF3PCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF3PCos_%i_%i", s3, s3, cent, ipt));

            hLamBarF2MCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF2MCos_%i_%i", s3, s3, cent, ipt));
            hLamBarF2Cos2_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF2Cos2_%i_%i", s3, s3, cent, ipt));
            hLamBarF3MCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF3MCos_%i_%i", s3, s3, cent, ipt));
            hLamBarF2PCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF2PCos_%i_%i", s3, s3, cent, ipt));
            hLamBarF3PCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF3PCos_%i_%i", s3, s3, cent, ipt));

            hLambdaPV2M_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaPV2M_%i_%i", s3, s3, cent, ipt));
            hLambdaMV2P_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaMV2P_%i_%i", s3, s3, cent, ipt));
            hLambdaFV2M_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaFV2M_%i_%i", s3, s3, cent, ipt));
            hLambdaFV2P_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaFV2P_%i_%i", s3, s3, cent, ipt));
            hLamBarPV2M_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarPV2M_%i_%i", s3, s3, cent, ipt));
            hLamBarMV2P_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarMV2P_%i_%i", s3, s3, cent, ipt));
            hLamBarFV2M_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarFV2M_%i_%i", s3, s3, cent, ipt));
            hLamBarFV2P_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarFV2P_%i_%i", s3, s3, cent, ipt));

            hLambdaPV3M_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaPV3M_%i_%i", s3, s3, cent, ipt));
            hLambdaMV3P_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaMV3P_%i_%i", s3, s3, cent, ipt));
            hLambdaFV3M_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaFV3M_%i_%i", s3, s3, cent, ipt));
            hLambdaFV3P_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaFV3P_%i_%i", s3, s3, cent, ipt));
            hLamBarPV3M_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarPV3M_%i_%i", s3, s3, cent, ipt));
            hLamBarMV3P_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarMV3P_%i_%i", s3, s3, cent, ipt));
            hLamBarFV3M_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarFV3M_%i_%i", s3, s3, cent, ipt));
            hLamBarFV3P_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarFV3P_%i_%i", s3, s3, cent, ipt));

            if (bDebug) {
                cout << " --> cent = " << cent << " ipt = " << ipt << endl;
                cout << " --> hLambdaPV3M_ [cent][ipt] = " << hLambdaPV3M_ [cent][ipt]  << endl;
            }
        }
    }

    ////
    // EP Res
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

    // SP Res
    hReHFp2HFm2SP_   ->Divide(hCent2_);
    hReHFp2TrkMid2SP_->Divide(hCent2_);
    hReHFm2TrkMid2SP_->Divide(hCent2_);
    hImHFp2HFm2SP_   ->Divide(hCent2_);
    hImHFp2TrkMid2SP_->Divide(hCent2_);
    hImHFm2TrkMid2SP_->Divide(hCent2_);

    hReHFp3HFm3SP_   ->Divide(hCent2_);
    hReHFp3TrkMid3SP_->Divide(hCent2_);
    hReHFm3TrkMid3SP_->Divide(hCent2_);
    hImHFp3HFm3SP_   ->Divide(hCent2_);
    hImHFp3TrkMid3SP_->Divide(hCent2_);
    hImHFm3TrkMid3SP_->Divide(hCent2_);

    // cosTheta^2 combine P+M
    vector<vector<TH1D*>>   hLambdaCos2_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarCos2_(NCent, vector<TH1D*>(NpT));

    // Polar combine P+M
    vector<vector<TH1D*>>   hLambda2Cos_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar2Cos_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambda3Cos_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar3Cos_(NCent, vector<TH1D*>(NpT));

    // Polar SP combine P+M
    vector<vector<TH1D*>>   hLambda2CosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar2CosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambda3CosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar3CosSP_(NCent, vector<TH1D*>(NpT));

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // cos^2
            hLambdaCos2_[cent][ipt] = (TH1D*) hLambdaP2Cos2_[cent][ipt]->Clone(Form("hLambdaCos2_%i_%i", cent, ipt));
            hLamBarCos2_[cent][ipt] = (TH1D*) hLamBarP2Cos2_[cent][ipt]->Clone(Form("hLamBarCos2_%i_%i", cent, ipt));
            hLambdaCos2_[cent][ipt]->Add(hLambdaM2Cos2_[cent][ipt]);
            hLamBarCos2_[cent][ipt]->Add(hLamBarM2Cos2_[cent][ipt]);
            hLambdaCos2_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBarCos2_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaP2Cos2_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM2Cos2_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP2Cos2_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM2Cos2_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            // 2nd EP
            hLambda2Cos_ [cent][ipt] = (TH1D*) hLambdaP2Cos_[cent][ipt]->Clone(Form("hLambda2Cos_%i_%i", cent, ipt));
            hLamBar2Cos_ [cent][ipt] = (TH1D*) hLamBarP2Cos_[cent][ipt]->Clone(Form("hLamBar2Cos_%i_%i", cent, ipt));
            hLambda2Cos_[cent][ipt]->Add(hLambdaM2Cos_[cent][ipt]);
            hLamBar2Cos_[cent][ipt]->Add(hLamBarM2Cos_[cent][ipt]);
            hLambda2Cos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBar2Cos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaP2Cos_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM2Cos_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP2Cos_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM2Cos_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            // 2nd SP
            hLambda2CosSP_ [cent][ipt] = (TH1D*) hLambdaP2CosSP_[cent][ipt]->Clone(Form("hLambda2CosSP_%i_%i", cent, ipt));
            hLamBar2CosSP_ [cent][ipt] = (TH1D*) hLamBarP2CosSP_[cent][ipt]->Clone(Form("hLamBar2CosSP_%i_%i", cent, ipt));
            hLambda2CosSP_[cent][ipt]->Add(hLambdaM2CosSP_[cent][ipt]);
            hLamBar2CosSP_[cent][ipt]->Add(hLamBarM2CosSP_[cent][ipt]);
            hLambda2CosSP_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBar2CosSP_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaP2CosSP_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM2CosSP_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP2CosSP_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM2CosSP_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            // 3rd EP
            hLambda3Cos_ [cent][ipt] = (TH1D*) hLambdaP3Cos_[cent][ipt]->Clone(Form("hLambda3Cos_%i_%i", cent, ipt));
            hLamBar3Cos_ [cent][ipt] = (TH1D*) hLamBarP3Cos_[cent][ipt]->Clone(Form("hLamBar3Cos_%i_%i", cent, ipt));
            hLambda3Cos_[cent][ipt]->Add(hLambdaM3Cos_[cent][ipt]);
            hLamBar3Cos_[cent][ipt]->Add(hLamBarM3Cos_[cent][ipt]);
            hLambda3Cos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBar3Cos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaP3Cos_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM3Cos_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP3Cos_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM3Cos_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            // 3rd SP
            hLambda3CosSP_ [cent][ipt] = (TH1D*) hLambdaP3CosSP_[cent][ipt]->Clone(Form("hLambda3CosSP_%i_%i", cent, ipt));
            hLamBar3CosSP_ [cent][ipt] = (TH1D*) hLamBarP3CosSP_[cent][ipt]->Clone(Form("hLamBar3CosSP_%i_%i", cent, ipt));
            hLambda3CosSP_[cent][ipt]->Add(hLambdaM3CosSP_[cent][ipt]);
            hLamBar3CosSP_[cent][ipt]->Add(hLamBarM3CosSP_[cent][ipt]);
            hLambda3CosSP_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBar3CosSP_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaP3CosSP_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM3CosSP_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP3CosSP_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM3CosSP_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            // full eta
            hLambdaF2MCos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF2Cos2_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF3MCos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF2PCos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF3PCos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);

            hLamBarF2MCos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF2Cos2_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF3MCos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF2PCos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF3PCos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            // v2
            hLambdaPV2M_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaMV2P_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarPV2M_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarMV2P_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            hLambdaPV3M_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaMV3P_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarPV3M_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarMV3P_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            hLambdaFV2M_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaFV2P_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBarFV2M_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarFV2P_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaFV3M_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaFV3P_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBarFV3M_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarFV3P_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
        }
    }

    for ( int c = bpPb?6:0; c < NCent; c++ ) {
        double rHFp2HFm2 = hReHFp2HFm2_->GetBinContent(c);
        double iHFp2HFm2 = hImHFp2HFm2_->GetBinContent(c);
        double rHFp2Trk2 = hReHFp2TrkMid2_->GetBinContent(c);
        double iHFp2Trk2 = hImHFp2TrkMid2_->GetBinContent(c);
        double rHFm2Trk2 = hReHFm2TrkMid2_->GetBinContent(c);
        double iHFm2Trk2 = hImHFm2TrkMid2_->GetBinContent(c);

        std::complex cHFp2HFm2(rHFp2HFm2, iHFp2HFm2);
        std::complex cHFp2Trk2(rHFp2Trk2, iHFp2Trk2);
        std::complex cHFm2Trk2(rHFm2Trk2, iHFm2Trk2);

        double rHFp3HFm3 = hReHFp3HFm3_->GetBinContent(c);
        double iHFp3HFm3 = hImHFp3HFm3_->GetBinContent(c);
        double rHFp3Trk3 = hReHFp3TrkMid3_->GetBinContent(c);
        double iHFp3Trk3 = hImHFp3TrkMid3_->GetBinContent(c);
        double rHFm3Trk3 = hReHFm3TrkMid3_->GetBinContent(c);
        double iHFm3Trk3 = hImHFm3TrkMid3_->GetBinContent(c);

        std::complex cHFp3HFm3(rHFp3HFm3, iHFp3HFm3);
        std::complex cHFp3Trk3(rHFp3Trk3, iHFp3Trk3);
        std::complex cHFm3Trk3(rHFm3Trk3, iHFm3Trk3);

        std::complex resHFp2 = cHFp2HFm2 * cHFp2Trk2 / cHFm2Trk2;
        std::complex resHFm2 = std::conj(cHFp2HFm2) * cHFm2Trk2 / cHFp2Trk2;
        std::complex resHFp3 = cHFp3HFm3 * cHFp3Trk3 / cHFm3Trk3;
        std::complex resHFm3 = std::conj(cHFp3HFm3) * cHFm3Trk3 / cHFp3Trk3;

        if ( resHFp2.real() < 0 ) cout << "  --> cent = " << c << " resHFp2.real() = " << resHFp2.real() << endl;
        if ( resHFm2.real() < 0 ) cout << "  --> cent = " << c << " resHFm2.real() = " << resHFm2.real() << endl;
        if ( resHFp3.real() < 0 ) cout << "  --> cent = " << c << " resHFp3.real() = " << resHFp3.real() << endl;
        if ( resHFm3.real() < 0 ) cout << "  --> cent = " << c << " resHFm3.real() = " << resHFm3.real() << endl;

        hHF2pRes_->SetBinContent( c, resHFp2.real() );
        hHF2mRes_->SetBinContent( c, resHFm2.real() );
        hHF3pRes_->SetBinContent( c, resHFp3.real() );
        hHF3mRes_->SetBinContent( c, resHFm3.real() );
    }

    for ( int c = bpPb?6:0; c < NCent; c++ ) {
        double rHFp2HFm2 = hReHFp2HFm2SP_->GetBinContent(c);
        double iHFp2HFm2 = hImHFp2HFm2SP_->GetBinContent(c);
        double rHFp2Trk2 = hReHFp2TrkMid2SP_->GetBinContent(c);
        double iHFp2Trk2 = hImHFp2TrkMid2SP_->GetBinContent(c);
        double rHFm2Trk2 = hReHFm2TrkMid2SP_->GetBinContent(c);
        double iHFm2Trk2 = hImHFm2TrkMid2SP_->GetBinContent(c);

        std::complex cHFp2HFm2(rHFp2HFm2, iHFp2HFm2);
        std::complex cHFp2Trk2(rHFp2Trk2, iHFp2Trk2);
        std::complex cHFm2Trk2(rHFm2Trk2, iHFm2Trk2);

        double rHFp3HFm3 = hReHFp3HFm3SP_->GetBinContent(c);
        double iHFp3HFm3 = hImHFp3HFm3SP_->GetBinContent(c);
        double rHFp3Trk3 = hReHFp3TrkMid3SP_->GetBinContent(c);
        double iHFp3Trk3 = hImHFp3TrkMid3SP_->GetBinContent(c);
        double rHFm3Trk3 = hReHFm3TrkMid3SP_->GetBinContent(c);
        double iHFm3Trk3 = hImHFm3TrkMid3SP_->GetBinContent(c);

        std::complex cHFp3HFm3(rHFp3HFm3, iHFp3HFm3);
        std::complex cHFp3Trk3(rHFp3Trk3, iHFp3Trk3);
        std::complex cHFm3Trk3(rHFm3Trk3, iHFm3Trk3);

        std::complex resHFp2 = cHFp2HFm2 * cHFp2Trk2 / cHFm2Trk2;
        std::complex resHFm2 = std::conj(cHFp2HFm2) * cHFm2Trk2 / cHFp2Trk2;
        std::complex resHFp3 = cHFp3HFm3 * cHFp3Trk3 / cHFm3Trk3;
        std::complex resHFm3 = std::conj(cHFp3HFm3) * cHFm3Trk3 / cHFp3Trk3;

        if ( resHFp2.real() < 0 ) cout << "  --> cent = " << c << "  SP  resHFp2.real() = " << resHFp2.real() << endl;
        if ( resHFm2.real() < 0 ) cout << "  --> cent = " << c << "  SP  resHFm2.real() = " << resHFm2.real() << endl;
        if ( resHFp3.real() < 0 ) cout << "  --> cent = " << c << "  SP  resHFp3.real() = " << resHFp3.real() << endl;
        if ( resHFm3.real() < 0 ) cout << "  --> cent = " << c << "  SP  resHFm3.real() = " << resHFm3.real() << endl;

        hHF2pResSP_->SetBinContent( c, resHFp2.real() );
        hHF2mResSP_->SetBinContent( c, resHFm2.real() );
        hHF3pResSP_->SetBinContent( c, resHFp3.real() );
        hHF3mResSP_->SetBinContent( c, resHFm3.real() );
    }

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaP2Cos_  [cent][ipt]->Scale(1./sqrt(hHF2mRes_->GetBinContent(cent+1)));
            hLambdaP3Cos_  [cent][ipt]->Scale(1./sqrt(hHF3mRes_->GetBinContent(cent+1)));
            hLambdaP2CosSP_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaP3CosSP_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));

            hLambdaM2Cos_  [cent][ipt]->Scale(1./sqrt(hHF2pRes_->GetBinContent(cent+1)));
            hLambdaM3Cos_  [cent][ipt]->Scale(1./sqrt(hHF3pRes_->GetBinContent(cent+1)));
            hLambdaM2CosSP_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLambdaM3CosSP_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));

            hLambdaF2MCos_ [cent][ipt]->Scale(1./sqrt(hHF2mRes_->GetBinContent(cent+1)));
            hLambdaF3MCos_ [cent][ipt]->Scale(1./sqrt(hHF3mRes_->GetBinContent(cent+1)));
            hLambdaF2PCos_ [cent][ipt]->Scale(1./sqrt(hHF2pRes_->GetBinContent(cent+1)));
            hLambdaF3PCos_ [cent][ipt]->Scale(1./sqrt(hHF3pRes_->GetBinContent(cent+1)));

            hLambdaF2MCosSP_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaF3MCosSP_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaF2PCosSP_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLambdaF3PCosSP_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));

            hLamBarP2Cos_  [cent][ipt]->Scale(1./sqrt(hHF2mRes_->GetBinContent(cent+1)));
            hLamBarP3Cos_  [cent][ipt]->Scale(1./sqrt(hHF3mRes_->GetBinContent(cent+1)));
            hLamBarP2CosSP_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarP3CosSP_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));

            hLamBarM2Cos_  [cent][ipt]->Scale(1./sqrt(hHF2pRes_->GetBinContent(cent+1)));
            hLamBarM3Cos_  [cent][ipt]->Scale(1./sqrt(hHF3pRes_->GetBinContent(cent+1)));
            hLamBarM2CosSP_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarM3CosSP_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));

            hLamBarF2MCos_ [cent][ipt]->Scale(1./sqrt(hHF2mRes_->GetBinContent(cent+1)));
            hLamBarF3MCos_ [cent][ipt]->Scale(1./sqrt(hHF3mRes_->GetBinContent(cent+1)));
            hLamBarF2PCos_ [cent][ipt]->Scale(1./sqrt(hHF2pRes_->GetBinContent(cent+1)));
            hLamBarF3PCos_ [cent][ipt]->Scale(1./sqrt(hHF3pRes_->GetBinContent(cent+1)));

            hLamBarF2MCosSP_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarF3MCosSP_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarF2PCosSP_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarF3PCosSP_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
        }
    }


    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaPV2M_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaMV2P_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarPV2M_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarMV2P_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));

            hLambdaPV3M_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaMV3P_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLamBarPV3M_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarMV3P_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));

            hLambdaFV2M_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaFV2P_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarFV2M_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarFV2P_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));

            hLambdaFV3M_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaFV3P_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLamBarFV3M_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarFV3P_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
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

    hHF2pRes_->Write("hHF2pRes");
    hHF2mRes_->Write("hHF2mRes");
    hHF3pRes_->Write("hHF3pRes");
    hHF3mRes_->Write("hHF3mRes");

    hReHFp2HFm2SP_   ->Write("hReHFp2HFm2SP");
    hReHFp2TrkMid2SP_->Write("hReHFp2TrkMid2SP");
    hReHFm2TrkMid2SP_->Write("hReHFm2TrkMid2SP");
    hImHFp2HFm2SP_   ->Write("hImHFp2HFm2SP");
    hImHFp2TrkMid2SP_->Write("hImHFp2TrkMid2SP");
    hImHFm2TrkMid2SP_->Write("hImHFm2TrkMid2SP");

    hReHFp3HFm3SP_   ->Write("hReHFp3HFm3SP");
    hReHFp3TrkMid3SP_->Write("hReHFp3TrkMid3SP");
    hReHFm3TrkMid3SP_->Write("hReHFm3TrkMid3SP");
    hImHFp3HFm3SP_   ->Write("hImHFp3HFm3SP");
    hImHFp3TrkMid3SP_->Write("hImHFp3TrkMid3SP");
    hImHFm3TrkMid3SP_->Write("hImHFm3TrkMid3SP");

    hHF2pResSP_->Write("hHF2pResSP");
    hHF2mResSP_->Write("hHF2mResSP");
    hHF3pResSP_->Write("hHF3pResSP");
    hHF3mResSP_->Write("hHF3mResSP");

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        hHF2pO_    [cent]->Write(Form("hHF2pO_%i", cent));
        hHF2mO_    [cent]->Write(Form("hHF2mO_%i", cent));
        hHF2pF_    [cent]->Write(Form("hHF2pF_%i", cent));
        hHF2mF_    [cent]->Write(Form("hHF2mF_%i", cent));
        hHF3pO_    [cent]->Write(Form("hHF3pO_%i", cent));
        hHF3mO_    [cent]->Write(Form("hHF3mO_%i", cent));
        hHF3pF_    [cent]->Write(Form("hHF3pF_%i", cent));
        hHF3mF_    [cent]->Write(Form("hHF3mF_%i", cent));

        hNLambda_  [cent]->Write(Form("hNLambda_%i", cent));
        hLambdaEta_[cent]->Write(Form("hLambdaEta_%i", cent));

        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaMass_[cent][ipt]->Write(Form("hLambdaMass_%i_%i", cent, ipt));
            hLamBarMass_[cent][ipt]->Write(Form("hLamBarMass_%i_%i", cent, ipt));

            hLambdaP2Cos_ [cent][ipt]->Write(Form("hLambdaP2Cos_%i_%i",  cent, ipt));
            hLambdaP2Cos2_[cent][ipt]->Write(Form("hLambdaP2Cos2_%i_%i", cent, ipt));
            hLamBarP2Cos_ [cent][ipt]->Write(Form("hLamBarP2Cos_%i_%i",  cent, ipt));
            hLamBarP2Cos2_[cent][ipt]->Write(Form("hLamBarP2Cos2_%i_%i", cent, ipt));

            hLambdaM2Cos_ [cent][ipt]->Write(Form("hLambdaM2Cos_%i_%i",  cent, ipt));
            hLambdaM2Cos2_[cent][ipt]->Write(Form("hLambdaM2Cos2_%i_%i", cent, ipt));
            hLamBarM2Cos_ [cent][ipt]->Write(Form("hLamBarM2Cos_%i_%i",  cent, ipt));
            hLamBarM2Cos2_[cent][ipt]->Write(Form("hLamBarM2Cos2_%i_%i", cent, ipt));

            hLambdaP3Cos_ [cent][ipt]->Write(Form("hLambdaP3Cos_%i_%i", cent, ipt));
            hLamBarP3Cos_ [cent][ipt]->Write(Form("hLamBarP3Cos_%i_%i", cent, ipt));
            hLambdaM3Cos_ [cent][ipt]->Write(Form("hLambdaM3Cos_%i_%i", cent, ipt));
            hLamBarM3Cos_ [cent][ipt]->Write(Form("hLamBarM3Cos_%i_%i", cent, ipt));

            hLambdaP2CosSP_ [cent][ipt]->Write(Form("hLambdaP2CosSP_%i_%i", cent, ipt));
            hLambdaP3CosSP_ [cent][ipt]->Write(Form("hLambdaP3CosSP_%i_%i", cent, ipt));
            hLambdaM2CosSP_ [cent][ipt]->Write(Form("hLambdaM2CosSP_%i_%i", cent, ipt));
            hLambdaM3CosSP_ [cent][ipt]->Write(Form("hLambdaM3CosSP_%i_%i", cent, ipt));
            hLambdaF2MCosSP_[cent][ipt]->Write(Form("hLambdaF2MCosSP_%i_%i", cent, ipt));
            hLambdaF3MCosSP_[cent][ipt]->Write(Form("hLambdaF3MCosSP_%i_%i", cent, ipt));
            hLambdaF2PCosSP_[cent][ipt]->Write(Form("hLambdaF2PCosSP_%i_%i", cent, ipt));
            hLambdaF3PCosSP_[cent][ipt]->Write(Form("hLambdaF3PCosSP_%i_%i", cent, ipt));

            hLamBarP2CosSP_ [cent][ipt]->Write(Form("hLamBarP2CosSP_%i_%i", cent, ipt));
            hLamBarP3CosSP_ [cent][ipt]->Write(Form("hLamBarP3CosSP_%i_%i", cent, ipt));
            hLamBarM2CosSP_ [cent][ipt]->Write(Form("hLamBarM2CosSP_%i_%i", cent, ipt));
            hLamBarM3CosSP_ [cent][ipt]->Write(Form("hLamBarM3CosSP_%i_%i", cent, ipt));
            hLamBarF2MCosSP_[cent][ipt]->Write(Form("hLamBarF2MCosSP_%i_%i", cent, ipt));
            hLamBarF3MCosSP_[cent][ipt]->Write(Form("hLamBarF3MCosSP_%i_%i", cent, ipt));
            hLamBarF2PCosSP_[cent][ipt]->Write(Form("hLamBarF2PCosSP_%i_%i", cent, ipt));
            hLamBarF3PCosSP_[cent][ipt]->Write(Form("hLamBarF3PCosSP_%i_%i", cent, ipt));

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

            hLambdaCos2_[cent][ipt]->Write(Form("hLambdaCos2_%i_%i", cent, ipt));
            hLamBarCos2_[cent][ipt]->Write(Form("hLamBarCos2_%i_%i", cent, ipt));
            hLambda2Cos_[cent][ipt]->Write(Form("hLambda2Cos_%i_%i", cent, ipt));
            hLamBar2Cos_[cent][ipt]->Write(Form("hLamBar2Cos_%i_%i", cent, ipt));
            hLambda3Cos_[cent][ipt]->Write(Form("hLambda3Cos_%i_%i", cent, ipt));
            hLamBar3Cos_[cent][ipt]->Write(Form("hLamBar3Cos_%i_%i", cent, ipt));

            // v2
            hLambdaPV2M_[cent][ipt]->Write(Form("hLambdaPV2M_%i_%i", cent, ipt));
            hLambdaMV2P_[cent][ipt]->Write(Form("hLambdaMV2P_%i_%i", cent, ipt));
            hLamBarPV2M_[cent][ipt]->Write(Form("hLamBarPV2M_%i_%i", cent, ipt));
            hLamBarMV2P_[cent][ipt]->Write(Form("hLamBarMV2P_%i_%i", cent, ipt));

            hLambdaPV3M_[cent][ipt]->Write(Form("hLambdaPV3M_%i_%i", cent, ipt));
            hLambdaMV3P_[cent][ipt]->Write(Form("hLambdaMV3P_%i_%i", cent, ipt));
            hLamBarPV3M_[cent][ipt]->Write(Form("hLamBarPV3M_%i_%i", cent, ipt));
            hLamBarMV3P_[cent][ipt]->Write(Form("hLamBarMV3P_%i_%i", cent, ipt));

            hLambdaFV2M_[cent][ipt]->Write(Form("hLambdaFV2M_%i_%i", cent, ipt));
            hLambdaFV2P_[cent][ipt]->Write(Form("hLambdaFV2P_%i_%i", cent, ipt));
            hLamBarFV2M_[cent][ipt]->Write(Form("hLamBarFV2M_%i_%i", cent, ipt));
            hLamBarFV2P_[cent][ipt]->Write(Form("hLamBarFV2P_%i_%i", cent, ipt));

            hLambdaFV3M_[cent][ipt]->Write(Form("hLambdaFV3M_%i_%i", cent, ipt));
            hLambdaFV3P_[cent][ipt]->Write(Form("hLambdaFV3P_%i_%i", cent, ipt));
            hLamBarFV3M_[cent][ipt]->Write(Form("hLamBarFV3M_%i_%i", cent, ipt));
            hLamBarFV3P_[cent][ipt]->Write(Form("hLamBarFV3P_%i_%i", cent, ipt));
        }
    }
}
