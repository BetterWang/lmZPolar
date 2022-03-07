#include "label.h"
#include "const.h"
#include "complex"
#include "helper.h"

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
    std::cout << "   ---> " << ftxt[s1] << endl;

    TFile* f = new TFile(Form("%s/output_%i.root", ftxt[s1], s3));

    //
    TH1D*   hCent_ = (TH1D*) f->Get(Form("output_%i_%i/hCent", s3, s3));
    TH1D*   hCent2_= (TH1D*) hCent_->Rebin(NCent, "hCent2", CentBins);

    // EP 5% binning
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
    }

    // EP Resolution Parts
    // Re - real
    // Im - imag
    TH1D* hReHFp2HFm2_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2HFm2", s3, s3));
    TH1D* hReHFp2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2TrkMid2", s3, s3));
    TH1D* hReHFm2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm2TrkMid2", s3, s3));
    TH1D* hImHFp2HFm2_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2HFm2", s3, s3));
    TH1D* hImHFp2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2TrkMid2", s3, s3));
    TH1D* hImHFm2TrkMid2_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm2TrkMid2", s3, s3));

    TH1D* hReHFp3HFm3_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3HFm3", s3, s3));
    TH1D* hReHFp3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3TrkMid3", s3, s3));
    TH1D* hReHFm3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm3TrkMid3", s3, s3));
    TH1D* hImHFp3HFm3_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3HFm3", s3, s3));
    TH1D* hImHFp3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3TrkMid3", s3, s3));
    TH1D* hImHFm3TrkMid3_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm3TrkMid3", s3, s3));

    cout << "---> " << __LINE__ << endl;
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

    cout << "---> " << __LINE__ << endl;
    TH1D* hHF2pRes_ = new TH1D("hHF2pRes_", "hHF2pRes_", NCent, CentBins);
    TH1D* hHF2mRes_ = new TH1D("hHF2mRes_", "hHF2mRes_", NCent, CentBins);
    TH1D* hHF3pRes_ = new TH1D("hHF3pRes_", "hHF3pRes_", NCent, CentBins);
    TH1D* hHF3mRes_ = new TH1D("hHF3mRes_", "hHF3mRes_", NCent, CentBins);

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

    // SP Resolution Parts
    // Re - real
    // Im - imag
    TH1D* hReHFp2HFm2SP_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2HFm2SP", s3, s3));
    TH1D* hReHFp2TrkMid2SP_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp2TrkMid2SP", s3, s3));
    TH1D* hReHFm2TrkMid2SP_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm2TrkMid2SP", s3, s3));
    TH1D* hImHFp2HFm2SP_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2HFm2SP", s3, s3));
    TH1D* hImHFp2TrkMid2SP_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp2TrkMid2SP", s3, s3));
    TH1D* hImHFm2TrkMid2SP_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm2TrkMid2SP", s3, s3));

    TH1D* hReHFp3HFm3SP_    = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3HFm3SP", s3, s3));
    TH1D* hReHFp3TrkMid3SP_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFp3TrkMid3SP", s3, s3));
    TH1D* hReHFm3TrkMid3SP_ = (TH1D*) f->Get(Form("output_%i_%i/hReHFm3TrkMid3SP", s3, s3));
    TH1D* hImHFp3HFm3SP_    = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3HFm3SP", s3, s3));
    TH1D* hImHFp3TrkMid3SP_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFp3TrkMid3SP", s3, s3));
    TH1D* hImHFm3TrkMid3SP_ = (TH1D*) f->Get(Form("output_%i_%i/hImHFm3TrkMid3SP", s3, s3));

    ////
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

    TH1D* hHF2pResSP_ = new TH1D("hHF2pResSP_", "hHF2pResSP_", NCent, CentBins);
    TH1D* hHF2mResSP_ = new TH1D("hHF2mResSP_", "hHF2mResSP_", NCent, CentBins);
    TH1D* hHF3pResSP_ = new TH1D("hHF3pResSP_", "hHF3pResSP_", NCent, CentBins);
    TH1D* hHF3mResSP_ = new TH1D("hHF3mResSP_", "hHF3mResSP_", NCent, CentBins);

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

    // Lm mass binning
    // [cent][ipt]
    vector<vector<TH1D*>>   hLambdaMass_   (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMass_   (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMassP_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMassP_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMassM_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMassM_  (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaReMassP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarReMassP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaReMassM_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarReMassM_(NCent, vector<TH1D*>(NpT));

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaMass_[cent][ipt]  = (TH1D*) f->Get(Form("output_%i_%i/hLambdaMass_%i_%i", s3, s3, cent, ipt));
            hLamBarMass_[cent][ipt]  = (TH1D*) f->Get(Form("output_%i_%i/hLamBarMass_%i_%i", s3, s3, cent, ipt));
            hLambdaMassP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaMassP_%i_%i", s3, s3, cent, ipt));
            hLamBarMassP_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarMassP_%i_%i", s3, s3, cent, ipt));
            hLambdaMassM_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaMassM_%i_%i", s3, s3, cent, ipt));
            hLamBarMassM_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarMassM_%i_%i", s3, s3, cent, ipt));

            hLambdaReMass_[cent][ipt]  = RebinLmMass(Form("hLambdaReMass_%i_%i",  cent, ipt ), hLambdaMass_[cent][ipt]  );
            hLamBarReMass_[cent][ipt]  = RebinLmMass(Form("hLamBarReMass_%i_%i",  cent, ipt ), hLamBarMass_[cent][ipt]  );
            hLambdaReMassP_[cent][ipt] = RebinLmMass(Form("hLambdaReMassP_%i_%i", cent, ipt ), hLambdaMassP_[cent][ipt] );
            hLamBarReMassP_[cent][ipt] = RebinLmMass(Form("hLamBarReMassP_%i_%i", cent, ipt ), hLamBarMassP_[cent][ipt] );
            hLambdaReMassM_[cent][ipt] = RebinLmMass(Form("hLambdaReMassM_%i_%i", cent, ipt ), hLambdaMassM_[cent][ipt] );
            hLamBarReMassM_[cent][ipt] = RebinLmMass(Form("hLamBarReMassM_%i_%i", cent, ipt ), hLamBarMassM_[cent][ipt] );
        }
    }

    // v2 v3
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

    vector<vector<TH1D*>>   hLambdaPV2MReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMV2PReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV2MReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV2PReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarPV2MReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMV2PReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV2MReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV2PReMass_ (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaPV3MReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMV3PReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV3MReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV3PReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarPV3MReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMV3PReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV3MReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV3PReMass_ (NCent, vector<TH1D*>(NpT));

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
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

            hLambdaPV2MReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaPV2MReMass_%i_%i", cent, ipt), hLambdaPV2M_[cent][ipt]);
            hLambdaMV2PReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaMV2PReMass_%i_%i", cent, ipt), hLambdaMV2P_[cent][ipt]);
            hLambdaFV2MReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaFV2MReMass_%i_%i", cent, ipt), hLambdaFV2M_[cent][ipt]);
            hLambdaFV2PReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaFV2PReMass_%i_%i", cent, ipt), hLambdaFV2P_[cent][ipt]);
            hLamBarPV2MReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarPV2MReMass_%i_%i", cent, ipt), hLamBarPV2M_[cent][ipt]);
            hLamBarMV2PReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarMV2PReMass_%i_%i", cent, ipt), hLamBarMV2P_[cent][ipt]);
            hLamBarFV2MReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarFV2MReMass_%i_%i", cent, ipt), hLamBarFV2M_[cent][ipt]);
            hLamBarFV2PReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarFV2PReMass_%i_%i", cent, ipt), hLamBarFV2P_[cent][ipt]);

            hLambdaPV3MReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaPV3MReMass_%i_%i", cent, ipt), hLambdaPV3M_[cent][ipt]);
            hLambdaMV3PReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaMV3PReMass_%i_%i", cent, ipt), hLambdaMV3P_[cent][ipt]);
            hLambdaFV3MReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaFV3MReMass_%i_%i", cent, ipt), hLambdaFV3M_[cent][ipt]);
            hLambdaFV3PReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaFV3PReMass_%i_%i", cent, ipt), hLambdaFV3P_[cent][ipt]);
            hLamBarPV3MReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarPV3MReMass_%i_%i", cent, ipt), hLamBarPV3M_[cent][ipt]);
            hLamBarMV3PReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarMV3PReMass_%i_%i", cent, ipt), hLamBarMV3P_[cent][ipt]);
            hLamBarFV3MReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarFV3MReMass_%i_%i", cent, ipt), hLamBarFV3M_[cent][ipt]);
            hLamBarFV3PReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarFV3PReMass_%i_%i", cent, ipt), hLamBarFV3P_[cent][ipt]);
        }
    }

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // v2 and v3 before EP eff
            hLambdaPV2M_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaMV2P_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarPV2M_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarMV2P_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);
            hLambdaFV2M_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaFV2P_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBarFV2M_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarFV2P_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaPV3M_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaMV3P_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarPV3M_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarMV3P_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);
            hLambdaFV3M_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaFV3P_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBarFV3M_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarFV3P_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            //// v2 v3 ReMass
            hLambdaPV2MReMass_[cent][ipt]->Divide(hLambdaReMassP_[cent][ipt]);
            hLambdaMV2PReMass_[cent][ipt]->Divide(hLambdaReMassM_[cent][ipt]);
            hLambdaFV2MReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLambdaFV2PReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLamBarPV2MReMass_[cent][ipt]->Divide(hLamBarReMassP_[cent][ipt]);
            hLamBarMV2PReMass_[cent][ipt]->Divide(hLamBarReMassM_[cent][ipt]);
            hLamBarFV2MReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLamBarFV2PReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);

            hLambdaPV3MReMass_[cent][ipt]->Divide(hLambdaReMassP_[cent][ipt]);
            hLambdaMV3PReMass_[cent][ipt]->Divide(hLambdaReMassM_[cent][ipt]);
            hLambdaFV3MReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLambdaFV3PReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLamBarPV3MReMass_[cent][ipt]->Divide(hLamBarReMassP_[cent][ipt]);
            hLamBarMV3PReMass_[cent][ipt]->Divide(hLamBarReMassM_[cent][ipt]);
            hLamBarFV3MReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLamBarFV3PReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);

            // v2 v3 EP res
            hLambdaPV2M_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaMV2P_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarPV2M_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarMV2P_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLambdaFV2M_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaFV2P_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarFV2M_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarFV2P_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));

            hLambdaPV3M_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaMV3P_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLamBarPV3M_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarMV3P_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLambdaFV3M_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaFV3P_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLamBarFV3M_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarFV3P_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));

            //// v2 v3 ReMass
            hLambdaPV2MReMass_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaMV2PReMass_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLambdaFV2MReMass_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaFV2PReMass_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarPV2MReMass_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarMV2PReMass_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarFV2MReMass_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarFV2PReMass_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));

            hLambdaPV3MReMass_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaMV3PReMass_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLambdaFV3MReMass_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaFV3PReMass_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLamBarPV3MReMass_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarMV3PReMass_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLamBarFV3MReMass_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarFV3PReMass_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
        }
    }

    // Polar cosTheta^2 -- acceptance term
    vector<vector<TH1D*>>   hLambdaP2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2Cos2_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2Cos2_ (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaP2Cos2ReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2Cos2ReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2Cos2ReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2Cos2ReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2Cos2ReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2Cos2ReMass_ (NCent, vector<TH1D*>(NpT));

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaP2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaP2Cos2_%i_%i", s3, s3, cent, ipt));
            hLambdaM2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaM2Cos2_%i_%i", s3, s3, cent, ipt));
            hLamBarP2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarP2Cos2_%i_%i", s3, s3, cent, ipt));
            hLamBarM2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarM2Cos2_%i_%i", s3, s3, cent, ipt));
            hLambdaF2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF2Cos2_%i_%i", s3, s3, cent, ipt));
            hLamBarF2Cos2_[cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF2Cos2_%i_%i", s3, s3, cent, ipt));

            // ReMass
            hLambdaP2Cos2ReMass_[cent][ipt] = RebinLmMass(Form("hLambdaP2Cos2ReMass_%i_%i", cent, ipt), hLambdaP2Cos2_[cent][ipt] );
            hLambdaM2Cos2ReMass_[cent][ipt] = RebinLmMass(Form("hLambdaM2Cos2ReMass_%i_%i", cent, ipt), hLambdaM2Cos2_[cent][ipt] );
            hLamBarP2Cos2ReMass_[cent][ipt] = RebinLmMass(Form("hLamBarP2Cos2ReMass_%i_%i", cent, ipt), hLamBarP2Cos2_[cent][ipt] );
            hLamBarM2Cos2ReMass_[cent][ipt] = RebinLmMass(Form("hLamBarM2Cos2ReMass_%i_%i", cent, ipt), hLamBarM2Cos2_[cent][ipt] );
            hLambdaF2Cos2ReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF2Cos2ReMass_%i_%i", cent, ipt), hLambdaF2Cos2_[cent][ipt] );
            hLamBarF2Cos2ReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF2Cos2ReMass_%i_%i", cent, ipt), hLamBarF2Cos2_[cent][ipt] );

            //
            hLambdaP2Cos2_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM2Cos2_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP2Cos2_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM2Cos2_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);
            hLambdaF2Cos2_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBarF2Cos2_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaP2Cos2ReMass_[cent][ipt]->Divide(hLambdaReMassP_[cent][ipt]);
            hLambdaM2Cos2ReMass_[cent][ipt]->Divide(hLambdaReMassM_[cent][ipt]);
            hLamBarP2Cos2ReMass_[cent][ipt]->Divide(hLamBarReMassP_[cent][ipt]);
            hLamBarM2Cos2ReMass_[cent][ipt]->Divide(hLamBarReMassM_[cent][ipt]);
            hLambdaF2Cos2ReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLamBarF2Cos2ReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
        }
    }

    // Polar P(M)-side wrt HFm(HFp)
    vector<vector<TH1D*>>   hLambdaP2Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2Cos_  (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaP3Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM3Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP3Cos_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM3Cos_  (NCent, vector<TH1D*>(NpT));

    //// full eta
    vector<vector<TH1D*>>   hLambdaF2MCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3MCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2PCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3PCos_ (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLamBarF2MCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3MCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2PCos_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3PCos_ (NCent, vector<TH1D*>(NpT));

    // Polar SP
    vector<vector<TH1D*>>   hLambdaP2CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2CosSP_ (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaP3CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM3CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP3CosSP_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM3CosSP_ (NCent, vector<TH1D*>(NpT));

    //// full eta
    vector<vector<TH1D*>>   hLambdaF2MCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2PCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2MCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2PCosSP_(NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaF3MCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3PCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3MCosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3PCosSP_(NCent, vector<TH1D*>(NpT));

    //// ReMass
    // Polar P(M)-side wrt HFm(HFp)
    vector<vector<TH1D*>>   hLambdaP2CosReMass_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2CosReMass_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2CosReMass_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2CosReMass_  (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaP3CosReMass_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM3CosReMass_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP3CosReMass_  (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM3CosReMass_  (NCent, vector<TH1D*>(NpT));

    //// full eta
    vector<vector<TH1D*>>   hLambdaF2MCosReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3MCosReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2PCosReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3PCosReMass_ (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLamBarF2MCosReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3MCosReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2PCosReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3PCosReMass_ (NCent, vector<TH1D*>(NpT));

    // Polar SP
    vector<vector<TH1D*>>   hLambdaP2CosSPReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2CosSPReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2CosSPReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2CosSPReMass_ (NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaP3CosSPReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM3CosSPReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP3CosSPReMass_ (NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM3CosSPReMass_ (NCent, vector<TH1D*>(NpT));

    //// full eta
    vector<vector<TH1D*>>   hLambdaF2MCosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2PCosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2MCosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2PCosSPReMass_(NCent, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaF3MCosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3PCosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3MCosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3PCosSPReMass_(NCent, vector<TH1D*>(NpT));

    // Get
    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // Polar P(M)-side wrt HFm(HFp)
            hLambdaP2Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaP2Cos_%i_%i", s3, s3, cent, ipt));
            hLambdaM2Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaM2Cos_%i_%i", s3, s3, cent, ipt));
            hLamBarP2Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarP2Cos_%i_%i", s3, s3, cent, ipt));
            hLamBarM2Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarM2Cos_%i_%i", s3, s3, cent, ipt));

            hLambdaP3Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaP3Cos_%i_%i", s3, s3, cent, ipt));
            hLambdaM3Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaM3Cos_%i_%i", s3, s3, cent, ipt));
            hLamBarP3Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarP3Cos_%i_%i", s3, s3, cent, ipt));
            hLamBarM3Cos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarM3Cos_%i_%i", s3, s3, cent, ipt));

            // full eta
            hLambdaF2MCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF2MCos_%i_%i", s3, s3, cent, ipt));
            hLambdaF3MCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF3MCos_%i_%i", s3, s3, cent, ipt));
            hLambdaF2PCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF2PCos_%i_%i", s3, s3, cent, ipt));
            hLambdaF3PCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLambdaF3PCos_%i_%i", s3, s3, cent, ipt));

            hLamBarF2MCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF2MCos_%i_%i", s3, s3, cent, ipt));
            hLamBarF3MCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF3MCos_%i_%i", s3, s3, cent, ipt));
            hLamBarF2PCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF2PCos_%i_%i", s3, s3, cent, ipt));
            hLamBarF3PCos_ [cent][ipt] = (TH1D*) f->Get(Form("output_%i_%i/hLamBarF3PCos_%i_%i", s3, s3, cent, ipt));


            // Polar SP
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
        }
    }

    // EP SP Res
    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaP2Cos_  [cent][ipt]->Scale(1./sqrt(hHF2mRes_->GetBinContent(cent+1)));
            hLambdaM2Cos_  [cent][ipt]->Scale(1./sqrt(hHF2pRes_->GetBinContent(cent+1)));
            hLamBarP2Cos_  [cent][ipt]->Scale(1./sqrt(hHF2mRes_->GetBinContent(cent+1)));
            hLamBarM2Cos_  [cent][ipt]->Scale(1./sqrt(hHF2pRes_->GetBinContent(cent+1)));

            hLambdaP3Cos_  [cent][ipt]->Scale(1./sqrt(hHF3mRes_->GetBinContent(cent+1)));
            hLambdaM3Cos_  [cent][ipt]->Scale(1./sqrt(hHF3pRes_->GetBinContent(cent+1)));
            hLamBarP3Cos_  [cent][ipt]->Scale(1./sqrt(hHF3mRes_->GetBinContent(cent+1)));
            hLamBarM3Cos_  [cent][ipt]->Scale(1./sqrt(hHF3pRes_->GetBinContent(cent+1)));

            hLambdaF2MCos_ [cent][ipt]->Scale(1./sqrt(hHF2mRes_->GetBinContent(cent+1)));
            hLambdaF3MCos_ [cent][ipt]->Scale(1./sqrt(hHF3mRes_->GetBinContent(cent+1)));
            hLambdaF2PCos_ [cent][ipt]->Scale(1./sqrt(hHF2pRes_->GetBinContent(cent+1)));
            hLambdaF3PCos_ [cent][ipt]->Scale(1./sqrt(hHF3pRes_->GetBinContent(cent+1)));

            hLamBarF2MCos_ [cent][ipt]->Scale(1./sqrt(hHF2mRes_->GetBinContent(cent+1)));
            hLamBarF3MCos_ [cent][ipt]->Scale(1./sqrt(hHF3mRes_->GetBinContent(cent+1)));
            hLamBarF2PCos_ [cent][ipt]->Scale(1./sqrt(hHF2pRes_->GetBinContent(cent+1)));
            hLamBarF3PCos_ [cent][ipt]->Scale(1./sqrt(hHF3pRes_->GetBinContent(cent+1)));

            hLambdaP2CosSP_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaM2CosSP_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarP2CosSP_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarM2CosSP_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));

            hLambdaP3CosSP_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaM3CosSP_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLamBarP3CosSP_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarM3CosSP_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));

            hLambdaF2MCosSP_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLambdaF2PCosSP_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));
            hLamBarF2MCosSP_[cent][ipt]->Scale(1./sqrt(hHF2mResSP_->GetBinContent(cent+1)));
            hLamBarF2PCosSP_[cent][ipt]->Scale(1./sqrt(hHF2pResSP_->GetBinContent(cent+1)));

            hLambdaF3MCosSP_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLambdaF3PCosSP_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
            hLamBarF3MCosSP_[cent][ipt]->Scale(1./sqrt(hHF3mResSP_->GetBinContent(cent+1)));
            hLamBarF3PCosSP_[cent][ipt]->Scale(1./sqrt(hHF3pResSP_->GetBinContent(cent+1)));
        }
    }

    // ReMass
    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // Polar P(M)-side wrt HFm(HFp)
            hLambdaP2CosReMass_[cent][ipt] = RebinLmMass(Form("hLambdaP2CosReMass_%i_%i", cent, ipt), hLambdaP2Cos_[cent][ipt] );
            hLambdaM2CosReMass_[cent][ipt] = RebinLmMass(Form("hLambdaM2CosReMass_%i_%i", cent, ipt), hLambdaM2Cos_[cent][ipt] );
            hLamBarP2CosReMass_[cent][ipt] = RebinLmMass(Form("hLamBarP2CosReMass_%i_%i", cent, ipt), hLamBarP2Cos_[cent][ipt] );
            hLamBarM2CosReMass_[cent][ipt] = RebinLmMass(Form("hLamBarM2CosReMass_%i_%i", cent, ipt), hLamBarM2Cos_[cent][ipt] );

            hLambdaP3CosReMass_[cent][ipt] = RebinLmMass(Form("hLambdaP3CosReMass_%i_%i", cent, ipt), hLambdaP3Cos_[cent][ipt] );
            hLambdaM3CosReMass_[cent][ipt] = RebinLmMass(Form("hLambdaM3CosReMass_%i_%i", cent, ipt), hLambdaM3Cos_[cent][ipt] );
            hLamBarP3CosReMass_[cent][ipt] = RebinLmMass(Form("hLamBarP3CosReMass_%i_%i", cent, ipt), hLamBarP3Cos_[cent][ipt] );
            hLamBarM3CosReMass_[cent][ipt] = RebinLmMass(Form("hLamBarM3CosReMass_%i_%i", cent, ipt), hLamBarM3Cos_[cent][ipt] );

            // full eta
            hLambdaF2MCosReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF2MCosReMass_%i_%i", cent, ipt), hLambdaF2MCos_[cent][ipt] );
            hLambdaF3MCosReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF3MCosReMass_%i_%i", cent, ipt), hLambdaF3MCos_[cent][ipt] );
            hLambdaF2PCosReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF2PCosReMass_%i_%i", cent, ipt), hLambdaF2PCos_[cent][ipt] );
            hLambdaF3PCosReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF3PCosReMass_%i_%i", cent, ipt), hLambdaF3PCos_[cent][ipt] );

            hLamBarF2MCosReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF2MCosReMass_%i_%i", cent, ipt), hLamBarF2MCos_[cent][ipt] );
            hLamBarF3MCosReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF3MCosReMass_%i_%i", cent, ipt), hLamBarF3MCos_[cent][ipt] );
            hLamBarF2PCosReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF2PCosReMass_%i_%i", cent, ipt), hLamBarF2PCos_[cent][ipt] );
            hLamBarF3PCosReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF3PCosReMass_%i_%i", cent, ipt), hLamBarF3PCos_[cent][ipt] );

            // Polar SP
            hLambdaP2CosSPReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaP2CosSPReMass_%i_%i",  cent, ipt), hLambdaP2CosSP_ [cent][ipt] );
            hLambdaM2CosSPReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaM2CosSPReMass_%i_%i",  cent, ipt), hLambdaM2CosSP_ [cent][ipt] );
            hLamBarP2CosSPReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarP2CosSPReMass_%i_%i",  cent, ipt), hLamBarP2CosSP_ [cent][ipt] );
            hLamBarM2CosSPReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarM2CosSPReMass_%i_%i",  cent, ipt), hLamBarM2CosSP_ [cent][ipt] );

            hLambdaP3CosSPReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaP3CosSPReMass_%i_%i",  cent, ipt), hLambdaP3CosSP_ [cent][ipt] );
            hLambdaM3CosSPReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaM3CosSPReMass_%i_%i",  cent, ipt), hLambdaM3CosSP_ [cent][ipt] );
            hLamBarP3CosSPReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarP3CosSPReMass_%i_%i",  cent, ipt), hLamBarP3CosSP_ [cent][ipt] );
            hLamBarM3CosSPReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarM3CosSPReMass_%i_%i",  cent, ipt), hLamBarM3CosSP_ [cent][ipt] );

            //// full eta
            hLambdaF2MCosSPReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF2MCosSPReMass_%i_%i", cent, ipt), hLambdaF2MCosSP_[cent][ipt] );
            hLambdaF2PCosSPReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF2PCosSPReMass_%i_%i", cent, ipt), hLambdaF2PCosSP_[cent][ipt] );
            hLamBarF2MCosSPReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF2MCosSPReMass_%i_%i", cent, ipt), hLamBarF2MCosSP_[cent][ipt] );
            hLamBarF2PCosSPReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF2PCosSPReMass_%i_%i", cent, ipt), hLamBarF2PCosSP_[cent][ipt] );

            hLambdaF3MCosSPReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF3MCosSPReMass_%i_%i", cent, ipt), hLambdaF3MCosSP_[cent][ipt] );
            hLambdaF3PCosSPReMass_[cent][ipt] = RebinLmMass(Form("hLambdaF3PCosSPReMass_%i_%i", cent, ipt), hLambdaF3PCosSP_[cent][ipt] );
            hLamBarF3MCosSPReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF3MCosSPReMass_%i_%i", cent, ipt), hLamBarF3MCosSP_[cent][ipt] );
            hLamBarF3PCosSPReMass_[cent][ipt] = RebinLmMass(Form("hLamBarF3PCosSPReMass_%i_%i", cent, ipt), hLamBarF3PCosSP_[cent][ipt] );
        }
    }

    // Polar combine P+M
    vector<vector<TH1D*>>   hLambda2Cos_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar2Cos_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambda3Cos_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar3Cos_(NCent, vector<TH1D*>(NpT));

    // ReMass
    vector<vector<TH1D*>>   hLambda2CosReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar2CosReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambda3CosReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar3CosReMass_(NCent, vector<TH1D*>(NpT));

    // Polar SP combine P+M
    vector<vector<TH1D*>>   hLambda2CosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar2CosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambda3CosSP_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar3CosSP_(NCent, vector<TH1D*>(NpT));

    // ReMass
    vector<vector<TH1D*>>   hLambda2CosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar2CosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambda3CosSPReMass_(NCent, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBar3CosSPReMass_(NCent, vector<TH1D*>(NpT));

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // EP
            hLambda2Cos_[cent][ipt] = (TH1D*) hLambdaP2Cos_[cent][ipt]->Clone(Form("hLambda2Cos_%i_%i", cent, ipt));
            hLamBar2Cos_[cent][ipt] = (TH1D*) hLamBarP2Cos_[cent][ipt]->Clone(Form("hLamBar2Cos_%i_%i", cent, ipt));
            hLambda2Cos_[cent][ipt]->Add(hLambdaM2Cos_[cent][ipt]);
            hLamBar2Cos_[cent][ipt]->Add(hLamBarM2Cos_[cent][ipt]);

            hLambda3Cos_[cent][ipt] = (TH1D*) hLambdaP3Cos_[cent][ipt]->Clone(Form("hLambda3Cos_%i_%i", cent, ipt));
            hLamBar3Cos_[cent][ipt] = (TH1D*) hLamBarP3Cos_[cent][ipt]->Clone(Form("hLamBar3Cos_%i_%i", cent, ipt));
            hLambda3Cos_[cent][ipt]->Add(hLambdaM3Cos_[cent][ipt]);
            hLamBar3Cos_[cent][ipt]->Add(hLamBarM3Cos_[cent][ipt]);

            // SP
            hLambda2CosSP_[cent][ipt] = (TH1D*) hLambdaP2CosSP_[cent][ipt]->Clone(Form("hLambda2CosSP_%i_%i", cent, ipt));
            hLamBar2CosSP_[cent][ipt] = (TH1D*) hLamBarP2CosSP_[cent][ipt]->Clone(Form("hLamBar2CosSP_%i_%i", cent, ipt));
            hLambda2CosSP_[cent][ipt]->Add(hLambdaM2CosSP_[cent][ipt]);
            hLamBar2CosSP_[cent][ipt]->Add(hLamBarM2CosSP_[cent][ipt]);

            hLambda3CosSP_[cent][ipt] = (TH1D*) hLambdaP3CosSP_[cent][ipt]->Clone(Form("hLambda3CosSP_%i_%i", cent, ipt));
            hLamBar3CosSP_[cent][ipt] = (TH1D*) hLamBarP3CosSP_[cent][ipt]->Clone(Form("hLamBar3CosSP_%i_%i", cent, ipt));
            hLambda3CosSP_[cent][ipt]->Add(hLambdaM3CosSP_[cent][ipt]);
            hLamBar3CosSP_[cent][ipt]->Add(hLamBarM3CosSP_[cent][ipt]);
        }
    }

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // ReMass
            hLambda2CosReMass_[cent][ipt] = RebinLmMass(Form("hLambda2CosReMass_%i_%i", cent, ipt), hLambda2Cos_[cent][ipt] );
            hLamBar2CosReMass_[cent][ipt] = RebinLmMass(Form("hLamBar2CosReMass_%i_%i", cent, ipt), hLamBar2Cos_[cent][ipt] );
            hLambda3CosReMass_[cent][ipt] = RebinLmMass(Form("hLambda3CosReMass_%i_%i", cent, ipt), hLambda3Cos_[cent][ipt] );
            hLamBar3CosReMass_[cent][ipt] = RebinLmMass(Form("hLamBar3CosReMass_%i_%i", cent, ipt), hLamBar3Cos_[cent][ipt] );

            hLambda2CosSPReMass_[cent][ipt] = RebinLmMass(Form("hLambda2CosSPReMass_%i_%i", cent, ipt), hLambda2CosSP_[cent][ipt] );
            hLamBar2CosSPReMass_[cent][ipt] = RebinLmMass(Form("hLamBar2CosSPReMass_%i_%i", cent, ipt), hLamBar2CosSP_[cent][ipt] );
            hLambda3CosSPReMass_[cent][ipt] = RebinLmMass(Form("hLambda3CosSPReMass_%i_%i", cent, ipt), hLambda3CosSP_[cent][ipt] );
            hLamBar3CosSPReMass_[cent][ipt] = RebinLmMass(Form("hLamBar3CosSPReMass_%i_%i", cent, ipt), hLamBar3CosSP_[cent][ipt] );
        }
    }

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // Norm
            hLambda2Cos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBar2Cos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLambda3Cos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBar3Cos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaP2Cos_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM2Cos_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP2Cos_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM2Cos_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            hLambdaP3Cos_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM3Cos_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP3Cos_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM3Cos_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            // full eta
            hLambdaF2MCos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF3MCos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF2PCos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF3PCos_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);

            hLamBarF2MCos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF3MCos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF2PCos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF3PCos_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            // Norm SP
            hLambda2CosSP_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBar2CosSP_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLambda3CosSP_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLamBar3CosSP_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            hLambdaP2CosSP_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM2CosSP_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP2CosSP_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM2CosSP_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            hLambdaP3CosSP_[cent][ipt]->Divide(hLambdaMassP_[cent][ipt]);
            hLambdaM3CosSP_[cent][ipt]->Divide(hLambdaMassM_[cent][ipt]);
            hLamBarP3CosSP_[cent][ipt]->Divide(hLamBarMassP_[cent][ipt]);
            hLamBarM3CosSP_[cent][ipt]->Divide(hLamBarMassM_[cent][ipt]);

            // full eta
            hLambdaF2MCosSP_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF3MCosSP_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF2PCosSP_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);
            hLambdaF3PCosSP_[cent][ipt]->Divide(hLambdaMass_[cent][ipt]);

            hLamBarF2MCosSP_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF3MCosSP_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF2PCosSP_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);
            hLamBarF3PCosSP_[cent][ipt]->Divide(hLamBarMass_[cent][ipt]);

            // ReMass Norm
            hLambda2CosReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLamBar2CosReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLambda3CosReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLamBar3CosReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);

            hLambdaP2CosReMass_[cent][ipt]->Divide(hLambdaReMassP_[cent][ipt]);
            hLambdaM2CosReMass_[cent][ipt]->Divide(hLambdaReMassM_[cent][ipt]);
            hLamBarP2CosReMass_[cent][ipt]->Divide(hLamBarReMassP_[cent][ipt]);
            hLamBarM2CosReMass_[cent][ipt]->Divide(hLamBarReMassM_[cent][ipt]);

            hLambdaP3CosReMass_[cent][ipt]->Divide(hLambdaReMassP_[cent][ipt]);
            hLambdaM3CosReMass_[cent][ipt]->Divide(hLambdaReMassM_[cent][ipt]);
            hLamBarP3CosReMass_[cent][ipt]->Divide(hLamBarReMassP_[cent][ipt]);
            hLamBarM3CosReMass_[cent][ipt]->Divide(hLamBarReMassM_[cent][ipt]);

            // full eta
            hLambdaF2MCosReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLambdaF3MCosReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLambdaF2PCosReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLambdaF3PCosReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);

            hLamBarF2MCosReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLamBarF3MCosReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLamBarF2PCosReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLamBarF3PCosReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);

            // ReMass Norm SP
            hLambda2CosSPReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLamBar2CosSPReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLambda3CosSPReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLamBar3CosSPReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);

            hLambdaP2CosSPReMass_[cent][ipt]->Divide(hLambdaReMassP_[cent][ipt]);
            hLambdaM2CosSPReMass_[cent][ipt]->Divide(hLambdaReMassM_[cent][ipt]);
            hLamBarP2CosSPReMass_[cent][ipt]->Divide(hLamBarReMassP_[cent][ipt]);
            hLamBarM2CosSPReMass_[cent][ipt]->Divide(hLamBarReMassM_[cent][ipt]);

            hLambdaP3CosSPReMass_[cent][ipt]->Divide(hLambdaReMassP_[cent][ipt]);
            hLambdaM3CosSPReMass_[cent][ipt]->Divide(hLambdaReMassM_[cent][ipt]);
            hLamBarP3CosSPReMass_[cent][ipt]->Divide(hLamBarReMassP_[cent][ipt]);
            hLamBarM3CosSPReMass_[cent][ipt]->Divide(hLamBarReMassM_[cent][ipt]);

            // full eta
            hLambdaF2MCosSPReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLambdaF3MCosSPReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLambdaF2PCosSPReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);
            hLambdaF3PCosSPReMass_[cent][ipt]->Divide(hLambdaReMass_[cent][ipt]);

            hLamBarF2MCosSPReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLamBarF3MCosSPReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLamBarF2PCosSPReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
            hLamBarF3PCosSPReMass_[cent][ipt]->Divide(hLamBarReMass_[cent][ipt]);
        }
    }

    // ReCent
    vector<vector<TH1D*>>   hLambdaPV2MReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMV2PReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV2MReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV2PReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarPV2MReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMV2PReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV2MReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV2PReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaPV3MReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaMV3PReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV3MReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaFV3PReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarPV3MReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarMV3PReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV3MReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarFV3PReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    // Polar cosTheta^2 -- acceptance term
    vector<vector<TH1D*>>   hLambdaP2Cos2ReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2Cos2ReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2Cos2ReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2Cos2ReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2Cos2ReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2Cos2ReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    // Polar P(M)-side wrt HFm(HFp)
    vector<vector<TH1D*>>   hLambdaP2CosReMassReCent_   (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2CosReMassReCent_   (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2CosReMassReCent_   (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2CosReMassReCent_   (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaP3CosReMassReCent_   (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM3CosReMassReCent_   (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP3CosReMassReCent_   (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM3CosReMassReCent_   (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    //// full eta
    vector<vector<TH1D*>>   hLambdaF2MCosReMassReCent_  (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3MCosReMassReCent_  (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2PCosReMassReCent_  (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3PCosReMassReCent_  (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLamBarF2MCosReMassReCent_  (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3MCosReMassReCent_  (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2PCosReMassReCent_  (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3PCosReMassReCent_  (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    // Polar SP
    vector<vector<TH1D*>>   hLambdaP2CosSPReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM2CosSPReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP2CosSPReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM2CosSPReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaP3CosSPReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaM3CosSPReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarP3CosSPReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarM3CosSPReMassReCent_ (bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    //// full eta
    vector<vector<TH1D*>>   hLambdaF2MCosSPReMassReCent_(bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF2PCosSPReMassReCent_(bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2MCosSPReMassReCent_(bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF2PCosSPReMassReCent_(bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    vector<vector<TH1D*>>   hLambdaF3MCosSPReMassReCent_(bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLambdaF3PCosSPReMassReCent_(bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3MCosSPReMassReCent_(bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));
    vector<vector<TH1D*>>   hLamBarF3PCosSPReMassReCent_(bpPb?NCent:NCentPbPb2018Rebin, vector<TH1D*>(NpT));

    if ( bpPb ) {
        for ( int cent = 6; cent < NCent; cent++ ) {
            for ( int ipt = 0; ipt < NpT; ipt++ ) {
                // v2 v3
                hLambdaPV2MReMassReCent_[cent][ipt] = hLambdaPV2MReMass_[cent][ipt];
                hLambdaMV2PReMassReCent_[cent][ipt] = hLambdaMV2PReMass_[cent][ipt];
                hLambdaFV2MReMassReCent_[cent][ipt] = hLambdaFV2MReMass_[cent][ipt];
                hLambdaFV2PReMassReCent_[cent][ipt] = hLambdaFV2PReMass_[cent][ipt];
                hLamBarPV2MReMassReCent_[cent][ipt] = hLamBarPV2MReMass_[cent][ipt];
                hLamBarMV2PReMassReCent_[cent][ipt] = hLamBarMV2PReMass_[cent][ipt];
                hLamBarFV2MReMassReCent_[cent][ipt] = hLamBarFV2MReMass_[cent][ipt];
                hLamBarFV2PReMassReCent_[cent][ipt] = hLamBarFV2PReMass_[cent][ipt];

                hLambdaPV3MReMassReCent_[cent][ipt] = hLambdaPV3MReMass_[cent][ipt];
                hLambdaMV3PReMassReCent_[cent][ipt] = hLambdaMV3PReMass_[cent][ipt];
                hLambdaFV3MReMassReCent_[cent][ipt] = hLambdaFV3MReMass_[cent][ipt];
                hLambdaFV3PReMassReCent_[cent][ipt] = hLambdaFV3PReMass_[cent][ipt];
                hLamBarPV3MReMassReCent_[cent][ipt] = hLamBarPV3MReMass_[cent][ipt];
                hLamBarMV3PReMassReCent_[cent][ipt] = hLamBarMV3PReMass_[cent][ipt];
                hLamBarFV3MReMassReCent_[cent][ipt] = hLamBarFV3MReMass_[cent][ipt];
                hLamBarFV3PReMassReCent_[cent][ipt] = hLamBarFV3PReMass_[cent][ipt];

                // Polar cosTheta^2 -- acceptance term
                hLambdaP2Cos2ReMassReCent_[cent][ipt] = hLambdaP2Cos2ReMass_[cent][ipt];
                hLambdaM2Cos2ReMassReCent_[cent][ipt] = hLambdaM2Cos2ReMass_[cent][ipt];
                hLamBarP2Cos2ReMassReCent_[cent][ipt] = hLamBarP2Cos2ReMass_[cent][ipt];
                hLamBarM2Cos2ReMassReCent_[cent][ipt] = hLamBarM2Cos2ReMass_[cent][ipt];
                hLambdaF2Cos2ReMassReCent_[cent][ipt] = hLambdaF2Cos2ReMass_[cent][ipt];
                hLamBarF2Cos2ReMassReCent_[cent][ipt] = hLamBarF2Cos2ReMass_[cent][ipt];

                // Polar P(M)-side wrt HFm(HFp)
                hLambdaP2CosReMassReCent_   [cent][ipt] = hLambdaP2CosReMass_   [cent][ipt];
                hLambdaM2CosReMassReCent_   [cent][ipt] = hLambdaM2CosReMass_   [cent][ipt];
                hLamBarP2CosReMassReCent_   [cent][ipt] = hLamBarP2CosReMass_   [cent][ipt];
                hLamBarM2CosReMassReCent_   [cent][ipt] = hLamBarM2CosReMass_   [cent][ipt];

                hLambdaP3CosReMassReCent_   [cent][ipt] = hLambdaP3CosReMass_   [cent][ipt];
                hLambdaM3CosReMassReCent_   [cent][ipt] = hLambdaM3CosReMass_   [cent][ipt];
                hLamBarP3CosReMassReCent_   [cent][ipt] = hLamBarP3CosReMass_   [cent][ipt];
                hLamBarM3CosReMassReCent_   [cent][ipt] = hLamBarM3CosReMass_   [cent][ipt];

                //// full eta[cent][ipt] ;
                hLambdaF2MCosReMassReCent_  [cent][ipt] = hLambdaF2MCosReMass_  [cent][ipt];
                hLambdaF3MCosReMassReCent_  [cent][ipt] = hLambdaF3MCosReMass_  [cent][ipt];
                hLambdaF2PCosReMassReCent_  [cent][ipt] = hLambdaF2PCosReMass_  [cent][ipt];
                hLambdaF3PCosReMassReCent_  [cent][ipt] = hLambdaF3PCosReMass_  [cent][ipt];

                hLamBarF2MCosReMassReCent_  [cent][ipt] = hLamBarF2MCosReMass_  [cent][ipt];
                hLamBarF3MCosReMassReCent_  [cent][ipt] = hLamBarF3MCosReMass_  [cent][ipt];
                hLamBarF2PCosReMassReCent_  [cent][ipt] = hLamBarF2PCosReMass_  [cent][ipt];
                hLamBarF3PCosReMassReCent_  [cent][ipt] = hLamBarF3PCosReMass_  [cent][ipt];

                // Polar SP[cent][ipt] ;
                hLambdaP2CosSPReMassReCent_ [cent][ipt] = hLambdaP2CosSPReMass_ [cent][ipt];
                hLambdaM2CosSPReMassReCent_ [cent][ipt] = hLambdaM2CosSPReMass_ [cent][ipt];
                hLamBarP2CosSPReMassReCent_ [cent][ipt] = hLamBarP2CosSPReMass_ [cent][ipt];
                hLamBarM2CosSPReMassReCent_ [cent][ipt] = hLamBarM2CosSPReMass_ [cent][ipt];

                hLambdaP3CosSPReMassReCent_ [cent][ipt] = hLambdaP3CosSPReMass_ [cent][ipt];
                hLambdaM3CosSPReMassReCent_ [cent][ipt] = hLambdaM3CosSPReMass_ [cent][ipt];
                hLamBarP3CosSPReMassReCent_ [cent][ipt] = hLamBarP3CosSPReMass_ [cent][ipt];
                hLamBarM3CosSPReMassReCent_ [cent][ipt] = hLamBarM3CosSPReMass_ [cent][ipt];

                //// full eta[cent][ipt] ;
                hLambdaF2MCosSPReMassReCent_[cent][ipt] = hLambdaF2MCosSPReMass_[cent][ipt];
                hLambdaF2PCosSPReMassReCent_[cent][ipt] = hLambdaF2PCosSPReMass_[cent][ipt];
                hLamBarF2MCosSPReMassReCent_[cent][ipt] = hLamBarF2MCosSPReMass_[cent][ipt];
                hLamBarF2PCosSPReMassReCent_[cent][ipt] = hLamBarF2PCosSPReMass_[cent][ipt];

                hLambdaF3MCosSPReMassReCent_[cent][ipt] = hLambdaF3MCosSPReMass_[cent][ipt];
                hLambdaF3PCosSPReMassReCent_[cent][ipt] = hLambdaF3PCosSPReMass_[cent][ipt];
                hLamBarF3MCosSPReMassReCent_[cent][ipt] = hLamBarF3MCosSPReMass_[cent][ipt];
                hLamBarF3PCosSPReMassReCent_[cent][ipt] = hLamBarF3PCosSPReMass_[cent][ipt];
            }
        }
    } else {
        auto CentMerge = [&](auto &from, auto &to) {
            for ( int cent = 0; cent < NCentPbPb2018Rebin; cent++ ) {
                for ( int ipt = 0; ipt < NpT; ipt++ ) {
                    to[cent][ipt] = new TH1D((TString(from[cent][ipt]->GetName())+"_ReCent").Data(), "", 5, Lm_mass_ranges);
                    for ( int i = 0; i < 5; i++ ) {
                        double sum[5]    = {};
                        double sum_w2[5] = {};
                        for ( int c : vectReCent[cent] ) {
                            double sigma2 = from[c][ipt]->GetBinError(i+1) * from[c][ipt]->GetBinError(i+1);
                            double cont  = from[c][ipt]->GetBinContent(i+1);

                            sum[i] += cont/sigma2;
                            sum_w2[i] += 1./sigma2;
                        }
                        if ( TMath::IsNaN(sum_w2[i]) ) {
                            to[cent][ipt]->SetBinContent(i+1, -999.);
                            to[cent][ipt]->SetBinError  (i+1, 9999.);
                        } else {
                            to[cent][ipt]->SetBinContent(i+1, sum[i]/sum_w2[i]);
                            to[cent][ipt]->SetBinError  (i+1, 1./sqrt(sum_w2[i]));
                        }
                    }
                }
            }
        };
        // v2 v3
        CentMerge(hLambdaPV2MReMass_, hLambdaPV2MReMassReCent_);
        CentMerge(hLambdaMV2PReMass_, hLambdaMV2PReMassReCent_);
        CentMerge(hLambdaFV2MReMass_, hLambdaFV2MReMassReCent_);
        CentMerge(hLambdaFV2PReMass_, hLambdaFV2PReMassReCent_);
        CentMerge(hLamBarPV2MReMass_, hLamBarPV2MReMassReCent_);
        CentMerge(hLamBarMV2PReMass_, hLamBarMV2PReMassReCent_);
        CentMerge(hLamBarFV2MReMass_, hLamBarFV2MReMassReCent_);
        CentMerge(hLamBarFV2PReMass_, hLamBarFV2PReMassReCent_);

        CentMerge(hLambdaPV3MReMass_, hLambdaPV3MReMassReCent_);
        CentMerge(hLambdaMV3PReMass_, hLambdaMV3PReMassReCent_);
        CentMerge(hLambdaFV3MReMass_, hLambdaFV3MReMassReCent_);
        CentMerge(hLambdaFV3PReMass_, hLambdaFV3PReMassReCent_);
        CentMerge(hLamBarPV3MReMass_, hLamBarPV3MReMassReCent_);
        CentMerge(hLamBarMV3PReMass_, hLamBarMV3PReMassReCent_);
        CentMerge(hLamBarFV3MReMass_, hLamBarFV3MReMassReCent_);
        CentMerge(hLamBarFV3PReMass_, hLamBarFV3PReMassReCent_);

        // Polar cosTheta^2 -- acceptance term
        CentMerge(hLambdaP2Cos2ReMass_, hLambdaP2Cos2ReMassReCent_);
        CentMerge(hLambdaM2Cos2ReMass_, hLambdaM2Cos2ReMassReCent_);
        CentMerge(hLamBarP2Cos2ReMass_, hLamBarP2Cos2ReMassReCent_);
        CentMerge(hLamBarM2Cos2ReMass_, hLamBarM2Cos2ReMassReCent_);
        CentMerge(hLambdaF2Cos2ReMass_, hLambdaF2Cos2ReMassReCent_);
        CentMerge(hLamBarF2Cos2ReMass_, hLamBarF2Cos2ReMassReCent_);

        // Polar P(M)-side wrt HFm(HFp)
        CentMerge(hLambdaP2CosReMass_    , hLambdaP2CosReMassReCent_   );
        CentMerge(hLambdaM2CosReMass_    , hLambdaM2CosReMassReCent_   );
        CentMerge(hLamBarP2CosReMass_    , hLamBarP2CosReMassReCent_   );
        CentMerge(hLamBarM2CosReMass_    , hLamBarM2CosReMassReCent_   );

        CentMerge(hLambdaP3CosReMass_    , hLambdaP3CosReMassReCent_   );
        CentMerge(hLambdaM3CosReMass_    , hLambdaM3CosReMassReCent_   );
        CentMerge(hLamBarP3CosReMass_    , hLamBarP3CosReMassReCent_   );
        CentMerge(hLamBarM3CosReMass_    , hLamBarM3CosReMassReCent_   );

        //// full eta[cent][ipt] ;
        CentMerge(hLambdaF2MCosReMass_   , hLambdaF2MCosReMassReCent_  );
        CentMerge(hLambdaF3MCosReMass_   , hLambdaF3MCosReMassReCent_  );
        CentMerge(hLambdaF2PCosReMass_   , hLambdaF2PCosReMassReCent_  );
        CentMerge(hLambdaF3PCosReMass_   , hLambdaF3PCosReMassReCent_  );

        CentMerge(hLamBarF2MCosReMass_   , hLamBarF2MCosReMassReCent_  );
        CentMerge(hLamBarF3MCosReMass_   , hLamBarF3MCosReMassReCent_  );
        CentMerge(hLamBarF2PCosReMass_   , hLamBarF2PCosReMassReCent_  );
        CentMerge(hLamBarF3PCosReMass_   , hLamBarF3PCosReMassReCent_  );

        // Polar SP[cent][ipt] ;
        CentMerge(hLambdaP2CosSPReMass_  , hLambdaP2CosSPReMassReCent_ );
        CentMerge(hLambdaM2CosSPReMass_  , hLambdaM2CosSPReMassReCent_ );
        CentMerge(hLamBarP2CosSPReMass_  , hLamBarP2CosSPReMassReCent_ );
        CentMerge(hLamBarM2CosSPReMass_  , hLamBarM2CosSPReMassReCent_ );

        CentMerge(hLambdaP3CosSPReMass_  , hLambdaP3CosSPReMassReCent_ );
        CentMerge(hLambdaM3CosSPReMass_  , hLambdaM3CosSPReMassReCent_ );
        CentMerge(hLamBarP3CosSPReMass_  , hLamBarP3CosSPReMassReCent_ );
        CentMerge(hLamBarM3CosSPReMass_  , hLamBarM3CosSPReMassReCent_ );

        //// full eta[cent][ipt] ;
        CentMerge(hLambdaF2MCosSPReMass_ , hLambdaF2MCosSPReMassReCent_);
        CentMerge(hLambdaF2PCosSPReMass_ , hLambdaF2PCosSPReMassReCent_);
        CentMerge(hLamBarF2MCosSPReMass_ , hLamBarF2MCosSPReMassReCent_);
        CentMerge(hLamBarF2PCosSPReMass_ , hLamBarF2PCosSPReMassReCent_);

        CentMerge(hLambdaF3MCosSPReMass_ , hLambdaF3MCosSPReMassReCent_);
        CentMerge(hLambdaF3PCosSPReMass_ , hLambdaF3PCosSPReMassReCent_);
        CentMerge(hLamBarF3MCosSPReMass_ , hLamBarF3MCosSPReMassReCent_);
        CentMerge(hLamBarF3PCosSPReMass_ , hLamBarF3PCosSPReMassReCent_);
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

    for ( int cent = (bpPb?6:0); cent < NCent; cent++ ) {
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
    }

    for ( int cent = (bpPb?6:0); cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaMass_[cent][ipt]->Write(Form("hLambdaMass_%i_%i", cent, ipt));
            hLamBarMass_[cent][ipt]->Write(Form("hLamBarMass_%i_%i", cent, ipt));
            hLambdaMassP_[cent][ipt]->Write(Form("hLambdaMassP_%i_%i", cent, ipt));
            hLamBarMassP_[cent][ipt]->Write(Form("hLamBarMassP_%i_%i", cent, ipt));
            hLambdaMassM_[cent][ipt]->Write(Form("hLambdaMassM_%i_%i", cent, ipt));
            hLamBarMassM_[cent][ipt]->Write(Form("hLamBarMassM_%i_%i", cent, ipt));

            hLambdaReMass_ [cent][ipt]->Write(Form("hLambdaReMass_%i_%i", cent, ipt));
            hLamBarReMass_ [cent][ipt]->Write(Form("hLamBarReMass_%i_%i", cent, ipt));
            hLambdaReMassP_[cent][ipt]->Write(Form("hLambdaReMassP_%i_%i", cent, ipt));
            hLamBarReMassP_[cent][ipt]->Write(Form("hLamBarReMassP_%i_%i", cent, ipt));
            hLambdaReMassM_[cent][ipt]->Write(Form("hLambdaReMassM_%i_%i", cent, ipt));
            hLamBarReMassM_[cent][ipt]->Write(Form("hLamBarReMassM_%i_%i", cent, ipt));

            // v2 v3
            hLambdaPV2MReMass_[cent][ipt]->Write(Form("hLambdaPV2MReMass_%i_%i", cent, ipt));
            hLambdaMV2PReMass_[cent][ipt]->Write(Form("hLambdaMV2PReMass_%i_%i", cent, ipt));
            hLamBarPV2MReMass_[cent][ipt]->Write(Form("hLamBarPV2MReMass_%i_%i", cent, ipt));
            hLamBarMV2PReMass_[cent][ipt]->Write(Form("hLamBarMV2PReMass_%i_%i", cent, ipt));
            hLambdaFV2MReMass_[cent][ipt]->Write(Form("hLambdaFV2MReMass_%i_%i", cent, ipt));
            hLambdaFV2PReMass_[cent][ipt]->Write(Form("hLambdaFV2PReMass_%i_%i", cent, ipt));
            hLamBarFV2MReMass_[cent][ipt]->Write(Form("hLamBarFV2MReMass_%i_%i", cent, ipt));
            hLamBarFV2PReMass_[cent][ipt]->Write(Form("hLamBarFV2PReMass_%i_%i", cent, ipt));

            hLambdaPV3MReMass_[cent][ipt]->Write(Form("hLambdaPV3MReMass_%i_%i", cent, ipt));
            hLambdaMV3PReMass_[cent][ipt]->Write(Form("hLambdaMV3PReMass_%i_%i", cent, ipt));
            hLamBarPV3MReMass_[cent][ipt]->Write(Form("hLamBarPV3MReMass_%i_%i", cent, ipt));
            hLamBarMV3PReMass_[cent][ipt]->Write(Form("hLamBarMV3PReMass_%i_%i", cent, ipt));
            hLambdaFV3MReMass_[cent][ipt]->Write(Form("hLambdaFV3MReMass_%i_%i", cent, ipt));
            hLambdaFV3PReMass_[cent][ipt]->Write(Form("hLambdaFV3PReMass_%i_%i", cent, ipt));
            hLamBarFV3MReMass_[cent][ipt]->Write(Form("hLamBarFV3MReMass_%i_%i", cent, ipt));
            hLamBarFV3PReMass_[cent][ipt]->Write(Form("hLamBarFV3PReMass_%i_%i", cent, ipt));

            // Polar cosTheta^2 -- acceptance term
            hLambdaP2Cos2ReMass_[cent][ipt]->Write(Form("hLambdaP2Cos2ReMass_%i_%i", cent, ipt));
            hLambdaM2Cos2ReMass_[cent][ipt]->Write(Form("hLambdaM2Cos2ReMass_%i_%i", cent, ipt));
            hLamBarP2Cos2ReMass_[cent][ipt]->Write(Form("hLamBarP2Cos2ReMass_%i_%i", cent, ipt));
            hLamBarM2Cos2ReMass_[cent][ipt]->Write(Form("hLamBarM2Cos2ReMass_%i_%i", cent, ipt));
            hLambdaF2Cos2ReMass_[cent][ipt]->Write(Form("hLambdaF2Cos2ReMass_%i_%i", cent, ipt));
            hLamBarF2Cos2ReMass_[cent][ipt]->Write(Form("hLamBarF2Cos2ReMass_%i_%i", cent, ipt));

            // Polar P(M)-side wrt HFm(HFp)
            hLambdaP2CosReMass_[cent][ipt]->Write(Form("hLambdaP2CosReMass_%i_%i", cent, ipt));
            hLambdaM2CosReMass_[cent][ipt]->Write(Form("hLambdaM2CosReMass_%i_%i", cent, ipt));
            hLamBarP2CosReMass_[cent][ipt]->Write(Form("hLamBarP2CosReMass_%i_%i", cent, ipt));
            hLamBarM2CosReMass_[cent][ipt]->Write(Form("hLamBarM2CosReMass_%i_%i", cent, ipt));

            hLambdaP3CosReMass_[cent][ipt]->Write(Form("hLambdaP3CosReMass_%i_%i", cent, ipt));
            hLambdaM3CosReMass_[cent][ipt]->Write(Form("hLambdaM3CosReMass_%i_%i", cent, ipt));
            hLamBarP3CosReMass_[cent][ipt]->Write(Form("hLamBarP3CosReMass_%i_%i", cent, ipt));
            hLamBarM3CosReMass_[cent][ipt]->Write(Form("hLamBarM3CosReMass_%i_%i", cent, ipt));

            hLambdaF2MCosReMass_[cent][ipt]->Write(Form("hLambdaF2MCosReMass_%i_%i", cent, ipt));
            hLambdaF3MCosReMass_[cent][ipt]->Write(Form("hLambdaF3MCosReMass_%i_%i", cent, ipt));
            hLambdaF2PCosReMass_[cent][ipt]->Write(Form("hLambdaF2PCosReMass_%i_%i", cent, ipt));
            hLambdaF3PCosReMass_[cent][ipt]->Write(Form("hLambdaF3PCosReMass_%i_%i", cent, ipt));

            hLamBarF2MCosReMass_[cent][ipt]->Write(Form("hLamBarF2MCosReMass_%i_%i", cent, ipt));
            hLamBarF3MCosReMass_[cent][ipt]->Write(Form("hLamBarF3MCosReMass_%i_%i", cent, ipt));
            hLamBarF2PCosReMass_[cent][ipt]->Write(Form("hLamBarF2PCosReMass_%i_%i", cent, ipt));
            hLamBarF3PCosReMass_[cent][ipt]->Write(Form("hLamBarF3PCosReMass_%i_%i", cent, ipt));

            // Polar SP
            hLambdaP2CosSPReMass_[cent][ipt]->Write(Form("hLambdaP2CosSPReMass_%i_%i", cent, ipt));
            hLambdaM2CosSPReMass_[cent][ipt]->Write(Form("hLambdaM2CosSPReMass_%i_%i", cent, ipt));
            hLamBarP2CosSPReMass_[cent][ipt]->Write(Form("hLamBarP2CosSPReMass_%i_%i", cent, ipt));
            hLamBarM2CosSPReMass_[cent][ipt]->Write(Form("hLamBarM2CosSPReMass_%i_%i", cent, ipt));

            hLambdaP3CosSPReMass_[cent][ipt]->Write(Form("hLambdaP3CosSPReMass_%i_%i", cent, ipt));
            hLambdaM3CosSPReMass_[cent][ipt]->Write(Form("hLambdaM3CosSPReMass_%i_%i", cent, ipt));
            hLamBarP3CosSPReMass_[cent][ipt]->Write(Form("hLamBarP3CosSPReMass_%i_%i", cent, ipt));
            hLamBarM3CosSPReMass_[cent][ipt]->Write(Form("hLamBarM3CosSPReMass_%i_%i", cent, ipt));

            hLambdaF2MCosSPReMass_[cent][ipt]->Write(Form("hLambdaF2MCosSPReMass_%i_%i", cent, ipt));
            hLambdaF3MCosSPReMass_[cent][ipt]->Write(Form("hLambdaF3MCosSPReMass_%i_%i", cent, ipt));
            hLambdaF2PCosSPReMass_[cent][ipt]->Write(Form("hLambdaF2PCosSPReMass_%i_%i", cent, ipt));
            hLambdaF3PCosSPReMass_[cent][ipt]->Write(Form("hLambdaF3PCosSPReMass_%i_%i", cent, ipt));

            hLamBarF2MCosSPReMass_[cent][ipt]->Write(Form("hLamBarF2MCosSPReMass_%i_%i", cent, ipt));
            hLamBarF3MCosSPReMass_[cent][ipt]->Write(Form("hLamBarF3MCosSPReMass_%i_%i", cent, ipt));
            hLamBarF2PCosSPReMass_[cent][ipt]->Write(Form("hLamBarF2PCosSPReMass_%i_%i", cent, ipt));
            hLamBarF3PCosSPReMass_[cent][ipt]->Write(Form("hLamBarF3PCosSPReMass_%i_%i", cent, ipt));

            // Polar combine P+M
            hLambda2CosReMass_[cent][ipt]->Write(Form("hLambda2CosReMass_%i_%i", cent, ipt));
            hLamBar2CosReMass_[cent][ipt]->Write(Form("hLamBar2CosReMass_%i_%i", cent, ipt));
            hLambda3CosReMass_[cent][ipt]->Write(Form("hLambda3CosReMass_%i_%i", cent, ipt));
            hLamBar3CosReMass_[cent][ipt]->Write(Form("hLamBar3CosReMass_%i_%i", cent, ipt));

            // SP
            hLambda2CosSPReMass_[cent][ipt]->Write(Form("hLambda2CosSPReMass_%i_%i", cent, ipt));
            hLamBar2CosSPReMass_[cent][ipt]->Write(Form("hLamBar2CosSPReMass_%i_%i", cent, ipt));
            hLambda3CosSPReMass_[cent][ipt]->Write(Form("hLambda3CosSPReMass_%i_%i", cent, ipt));
            hLamBar3CosSPReMass_[cent][ipt]->Write(Form("hLamBar3CosSPReMass_%i_%i", cent, ipt));
        }
    }

    for ( int cent = (bpPb?6:0); cent < (bpPb?NCent:NCentPbPb2018Rebin); cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // v2 v3
            hLambdaPV2MReMassReCent_[cent][ipt]->Write(Form("hLambdaPV2MReMassReCent_%i_%i", cent, ipt));
            hLambdaMV2PReMassReCent_[cent][ipt]->Write(Form("hLambdaMV2PReMassReCent_%i_%i", cent, ipt));
            hLambdaFV2MReMassReCent_[cent][ipt]->Write(Form("hLambdaFV2MReMassReCent_%i_%i", cent, ipt));
            hLambdaFV2PReMassReCent_[cent][ipt]->Write(Form("hLambdaFV2PReMassReCent_%i_%i", cent, ipt));
            hLamBarPV2MReMassReCent_[cent][ipt]->Write(Form("hLamBarPV2MReMassReCent_%i_%i", cent, ipt));
            hLamBarMV2PReMassReCent_[cent][ipt]->Write(Form("hLamBarMV2PReMassReCent_%i_%i", cent, ipt));
            hLamBarFV2MReMassReCent_[cent][ipt]->Write(Form("hLamBarFV2MReMassReCent_%i_%i", cent, ipt));
            hLamBarFV2PReMassReCent_[cent][ipt]->Write(Form("hLamBarFV2PReMassReCent_%i_%i", cent, ipt));

            hLambdaPV3MReMassReCent_[cent][ipt]->Write(Form("hLambdaPV3MReMassReCent_%i_%i", cent, ipt));
            hLambdaMV3PReMassReCent_[cent][ipt]->Write(Form("hLambdaMV3PReMassReCent_%i_%i", cent, ipt));
            hLambdaFV3MReMassReCent_[cent][ipt]->Write(Form("hLambdaFV3MReMassReCent_%i_%i", cent, ipt));
            hLambdaFV3PReMassReCent_[cent][ipt]->Write(Form("hLambdaFV3PReMassReCent_%i_%i", cent, ipt));
            hLamBarPV3MReMassReCent_[cent][ipt]->Write(Form("hLamBarPV3MReMassReCent_%i_%i", cent, ipt));
            hLamBarMV3PReMassReCent_[cent][ipt]->Write(Form("hLamBarMV3PReMassReCent_%i_%i", cent, ipt));
            hLamBarFV3MReMassReCent_[cent][ipt]->Write(Form("hLamBarFV3MReMassReCent_%i_%i", cent, ipt));
            hLamBarFV3PReMassReCent_[cent][ipt]->Write(Form("hLamBarFV3PReMassReCent_%i_%i", cent, ipt));

            // Polar cosTheta^2 -- acceptance term
            hLambdaP2Cos2ReMassReCent_[cent][ipt]->Write(Form("hLambdaP2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLambdaM2Cos2ReMassReCent_[cent][ipt]->Write(Form("hLambdaM2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLamBarP2Cos2ReMassReCent_[cent][ipt]->Write(Form("hLamBarP2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLamBarM2Cos2ReMassReCent_[cent][ipt]->Write(Form("hLamBarM2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLambdaF2Cos2ReMassReCent_[cent][ipt]->Write(Form("hLambdaF2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLamBarF2Cos2ReMassReCent_[cent][ipt]->Write(Form("hLamBarF2Cos2ReMassReCent_%i_%i", cent, ipt));

            // Polar P(M)-side wrt HFm(HFp)
            hLambdaP2CosReMassReCent_   [cent][ipt]->Write(Form("hLambdaP2CosReMassReCent_%i_%i", cent, ipt));
            hLambdaM2CosReMassReCent_   [cent][ipt]->Write(Form("hLambdaM2CosReMassReCent_%i_%i", cent, ipt));
            hLamBarP2CosReMassReCent_   [cent][ipt]->Write(Form("hLamBarP2CosReMassReCent_%i_%i", cent, ipt));
            hLamBarM2CosReMassReCent_   [cent][ipt]->Write(Form("hLamBarM2CosReMassReCent_%i_%i", cent, ipt));

            hLambdaP3CosReMassReCent_   [cent][ipt]->Write(Form("hLambdaP3CosReMassReCent_%i_%i", cent, ipt));
            hLambdaM3CosReMassReCent_   [cent][ipt]->Write(Form("hLambdaM3CosReMassReCent_%i_%i", cent, ipt));
            hLamBarP3CosReMassReCent_   [cent][ipt]->Write(Form("hLamBarP3CosReMassReCent_%i_%i", cent, ipt));
            hLamBarM3CosReMassReCent_   [cent][ipt]->Write(Form("hLamBarM3CosReMassReCent_%i_%i", cent, ipt));

            //// full eta
            hLambdaF2MCosReMassReCent_  [cent][ipt]->Write(Form("hLambdaF2MCosReMassReCent_%i_%i", cent, ipt));
            hLambdaF3MCosReMassReCent_  [cent][ipt]->Write(Form("hLambdaF3MCosReMassReCent_%i_%i", cent, ipt));
            hLambdaF2PCosReMassReCent_  [cent][ipt]->Write(Form("hLambdaF2PCosReMassReCent_%i_%i", cent, ipt));
            hLambdaF3PCosReMassReCent_  [cent][ipt]->Write(Form("hLambdaF3PCosReMassReCent_%i_%i", cent, ipt));

            hLamBarF2MCosReMassReCent_  [cent][ipt]->Write(Form("hLamBarF2MCosReMassReCent_%i_%i", cent, ipt));
            hLamBarF3MCosReMassReCent_  [cent][ipt]->Write(Form("hLamBarF3MCosReMassReCent_%i_%i", cent, ipt));
            hLamBarF2PCosReMassReCent_  [cent][ipt]->Write(Form("hLamBarF2PCosReMassReCent_%i_%i", cent, ipt));
            hLamBarF3PCosReMassReCent_  [cent][ipt]->Write(Form("hLamBarF3PCosReMassReCent_%i_%i", cent, ipt));

            // Polar SP
            hLambdaP2CosSPReMassReCent_ [cent][ipt]->Write(Form("hLambdaP2CosSPReMassReCent_%i_%i", cent, ipt));
            hLambdaM2CosSPReMassReCent_ [cent][ipt]->Write(Form("hLambdaM2CosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarP2CosSPReMassReCent_ [cent][ipt]->Write(Form("hLamBarP2CosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarM2CosSPReMassReCent_ [cent][ipt]->Write(Form("hLamBarM2CosSPReMassReCent_%i_%i", cent, ipt));

            hLambdaP3CosSPReMassReCent_ [cent][ipt]->Write(Form("hLambdaP3CosSPReMassReCent_%i_%i", cent, ipt));
            hLambdaM3CosSPReMassReCent_ [cent][ipt]->Write(Form("hLambdaM3CosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarP3CosSPReMassReCent_ [cent][ipt]->Write(Form("hLamBarP3CosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarM3CosSPReMassReCent_ [cent][ipt]->Write(Form("hLamBarM3CosSPReMassReCent_%i_%i", cent, ipt));

            //// full eta
            hLambdaF2MCosSPReMassReCent_[cent][ipt]->Write(Form("hLambdaF2MCosSPReMassReCent_%i_%i", cent, ipt));
            hLambdaF2PCosSPReMassReCent_[cent][ipt]->Write(Form("hLambdaF2PCosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarF2MCosSPReMassReCent_[cent][ipt]->Write(Form("hLamBarF2MCosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarF2PCosSPReMassReCent_[cent][ipt]->Write(Form("hLamBarF2PCosSPReMassReCent_%i_%i", cent, ipt));

            hLambdaF3MCosSPReMassReCent_[cent][ipt]->Write(Form("hLambdaF3MCosSPReMassReCent_%i_%i", cent, ipt));
            hLambdaF3PCosSPReMassReCent_[cent][ipt]->Write(Form("hLambdaF3PCosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarF3MCosSPReMassReCent_[cent][ipt]->Write(Form("hLamBarF3MCosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarF3PCosSPReMassReCent_[cent][ipt]->Write(Form("hLamBarF3PCosSPReMassReCent_%i_%i", cent, ipt));

        }
    }

}

