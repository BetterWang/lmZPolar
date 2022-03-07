#include "label.h"
#include "const.h"
#include "helper.h"
#include <TH1D.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TMath.h>

void genPlot(int s1 = 0)
{
    std::cout << " s1 = " << s1 << std::endl;
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

    TFile * f = new TFile(Form("%s/bGet_%i.root", ftxt[s1], 10));

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

    // Get
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

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaMass_[cent][ipt]  = (TH1D*) f->Get(Form("hLambdaMass_%i_%i",  cent, ipt));
            hLamBarMass_[cent][ipt]  = (TH1D*) f->Get(Form("hLamBarMass_%i_%i",  cent, ipt));
            hLambdaMassP_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaMassP_%i_%i", cent, ipt));
            hLamBarMassP_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarMassP_%i_%i", cent, ipt));
            hLambdaMassM_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaMassM_%i_%i", cent, ipt));
            hLamBarMassM_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarMassM_%i_%i", cent, ipt));

            hLambdaReMass_[cent][ipt]  = (TH1D*) f->Get(Form("hLambdaReMass_%i_%i",  cent, ipt));
            hLamBarReMass_[cent][ipt]  = (TH1D*) f->Get(Form("hLamBarReMass_%i_%i",  cent, ipt));
            hLambdaReMassP_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaReMassP_%i_%i", cent, ipt));
            hLamBarReMassP_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarReMassP_%i_%i", cent, ipt));
            hLambdaReMassM_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaReMassM_%i_%i", cent, ipt));
            hLamBarReMassM_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarReMassM_%i_%i", cent, ipt));
        }
    }

    for ( int cent = (bpPb?6:0); cent < (bpPb?NCent:NCentPbPb2018Rebin); cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            // v2 v3
            hLambdaPV2MReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaPV2MReMassReCent_%i_%i", cent, ipt));
            hLambdaMV2PReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaMV2PReMassReCent_%i_%i", cent, ipt));
            hLambdaFV2MReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaFV2MReMassReCent_%i_%i", cent, ipt));
            hLambdaFV2PReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaFV2PReMassReCent_%i_%i", cent, ipt));
            hLamBarPV2MReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarPV2MReMassReCent_%i_%i", cent, ipt));
            hLamBarMV2PReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarMV2PReMassReCent_%i_%i", cent, ipt));
            hLamBarFV2MReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarFV2MReMassReCent_%i_%i", cent, ipt));
            hLamBarFV2PReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarFV2PReMassReCent_%i_%i", cent, ipt));

            hLambdaPV3MReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaPV3MReMassReCent_%i_%i", cent, ipt));
            hLambdaMV3PReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaMV3PReMassReCent_%i_%i", cent, ipt));
            hLambdaFV3MReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaFV3MReMassReCent_%i_%i", cent, ipt));
            hLambdaFV3PReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaFV3PReMassReCent_%i_%i", cent, ipt));
            hLamBarPV3MReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarPV3MReMassReCent_%i_%i", cent, ipt));
            hLamBarMV3PReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarMV3PReMassReCent_%i_%i", cent, ipt));
            hLamBarFV3MReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarFV3MReMassReCent_%i_%i", cent, ipt));
            hLamBarFV3PReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarFV3PReMassReCent_%i_%i", cent, ipt));

            // Polar cosTheta^2 -- acceptance term
            hLambdaP2Cos2ReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaP2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLambdaM2Cos2ReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaM2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLamBarP2Cos2ReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarP2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLamBarM2Cos2ReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarM2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLambdaF2Cos2ReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaF2Cos2ReMassReCent_%i_%i", cent, ipt));
            hLamBarF2Cos2ReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarF2Cos2ReMassReCent_%i_%i", cent, ipt));

            // Polar P(M)-side wrt HFm(HFp)
            hLambdaP2CosReMassReCent_   [cent][ipt] = (TH1D*) f->Get(Form("hLambdaP2CosReMassReCent_%i_%i", cent, ipt));
            hLambdaM2CosReMassReCent_   [cent][ipt] = (TH1D*) f->Get(Form("hLambdaM2CosReMassReCent_%i_%i", cent, ipt));
            hLamBarP2CosReMassReCent_   [cent][ipt] = (TH1D*) f->Get(Form("hLamBarP2CosReMassReCent_%i_%i", cent, ipt));
            hLamBarM2CosReMassReCent_   [cent][ipt] = (TH1D*) f->Get(Form("hLamBarM2CosReMassReCent_%i_%i", cent, ipt));

            hLambdaP3CosReMassReCent_   [cent][ipt] = (TH1D*) f->Get(Form("hLambdaP3CosReMassReCent_%i_%i", cent, ipt));
            hLambdaM3CosReMassReCent_   [cent][ipt] = (TH1D*) f->Get(Form("hLambdaM3CosReMassReCent_%i_%i", cent, ipt));
            hLamBarP3CosReMassReCent_   [cent][ipt] = (TH1D*) f->Get(Form("hLamBarP3CosReMassReCent_%i_%i", cent, ipt));
            hLamBarM3CosReMassReCent_   [cent][ipt] = (TH1D*) f->Get(Form("hLamBarM3CosReMassReCent_%i_%i", cent, ipt));

            //// full eta
            hLambdaF2MCosReMassReCent_  [cent][ipt] = (TH1D*) f->Get(Form("hLambdaF2MCosReMassReCent_%i_%i", cent, ipt));
            hLambdaF3MCosReMassReCent_  [cent][ipt] = (TH1D*) f->Get(Form("hLambdaF3MCosReMassReCent_%i_%i", cent, ipt));
            hLambdaF2PCosReMassReCent_  [cent][ipt] = (TH1D*) f->Get(Form("hLambdaF2PCosReMassReCent_%i_%i", cent, ipt));
            hLambdaF3PCosReMassReCent_  [cent][ipt] = (TH1D*) f->Get(Form("hLambdaF3PCosReMassReCent_%i_%i", cent, ipt));

            hLamBarF2MCosReMassReCent_  [cent][ipt] = (TH1D*) f->Get(Form("hLamBarF2MCosReMassReCent_%i_%i", cent, ipt));
            hLamBarF3MCosReMassReCent_  [cent][ipt] = (TH1D*) f->Get(Form("hLamBarF3MCosReMassReCent_%i_%i", cent, ipt));
            hLamBarF2PCosReMassReCent_  [cent][ipt] = (TH1D*) f->Get(Form("hLamBarF2PCosReMassReCent_%i_%i", cent, ipt));
            hLamBarF3PCosReMassReCent_  [cent][ipt] = (TH1D*) f->Get(Form("hLamBarF3PCosReMassReCent_%i_%i", cent, ipt));

            // Polar SP
            hLambdaP2CosSPReMassReCent_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaP2CosSPReMassReCent_%i_%i", cent, ipt));
            hLambdaM2CosSPReMassReCent_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaM2CosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarP2CosSPReMassReCent_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarP2CosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarM2CosSPReMassReCent_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarM2CosSPReMassReCent_%i_%i", cent, ipt));

            hLambdaP3CosSPReMassReCent_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaP3CosSPReMassReCent_%i_%i", cent, ipt));
            hLambdaM3CosSPReMassReCent_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaM3CosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarP3CosSPReMassReCent_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarP3CosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarM3CosSPReMassReCent_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarM3CosSPReMassReCent_%i_%i", cent, ipt));

            //// full eta
            hLambdaF2MCosSPReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaF2MCosSPReMassReCent_%i_%i", cent, ipt));
            hLambdaF2PCosSPReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaF2PCosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarF2MCosSPReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarF2MCosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarF2PCosSPReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarF2PCosSPReMassReCent_%i_%i", cent, ipt));

            hLambdaF3MCosSPReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaF3MCosSPReMassReCent_%i_%i", cent, ipt));
            hLambdaF3PCosSPReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLambdaF3PCosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarF3MCosSPReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarF3MCosSPReMassReCent_%i_%i", cent, ipt));
            hLamBarF3PCosSPReMassReCent_[cent][ipt] = (TH1D*) f->Get(Form("hLamBarF3PCosSPReMassReCent_%i_%i", cent, ipt));

        }
    }

    // Get Vn
    GrData grLambdaPV2M(bpPb, hLambdaPV2MReMassReCent_, "grLambdaPV2M");
    GrData grLambdaMV2P(bpPb, hLambdaMV2PReMassReCent_, "grLambdaMV2P");
    GrData grLamBarPV2M(bpPb, hLamBarPV2MReMassReCent_, "grLamBarPV2M");
    GrData grLamBarMV2P(bpPb, hLamBarMV2PReMassReCent_, "grLamBarMV2P");
    GrData grLambdaFV2M(bpPb, hLambdaFV2MReMassReCent_, "grLambdaFV2M");
    GrData grLambdaFV2P(bpPb, hLambdaFV2PReMassReCent_, "grLambdaFV2P");
    GrData grLamBarFV2M(bpPb, hLamBarFV2MReMassReCent_, "grLamBarFV2M");
    GrData grLamBarFV2P(bpPb, hLamBarFV2PReMassReCent_, "grLamBarFV2P");

    GrData grLambdaPV3M(bpPb, hLambdaPV3MReMassReCent_, "grLambdaPV3M");
    GrData grLambdaMV3P(bpPb, hLambdaMV3PReMassReCent_, "grLambdaMV3P");
    GrData grLamBarPV3M(bpPb, hLamBarPV3MReMassReCent_, "grLamBarPV3M");
    GrData grLamBarMV3P(bpPb, hLamBarMV3PReMassReCent_, "grLamBarMV3P");
    GrData grLambdaFV3M(bpPb, hLambdaFV3MReMassReCent_, "grLambdaFV3M");
    GrData grLambdaFV3P(bpPb, hLambdaFV3PReMassReCent_, "grLambdaFV3P");
    GrData grLamBarFV3M(bpPb, hLamBarFV3MReMassReCent_, "grLamBarFV3M");
    GrData grLamBarFV3P(bpPb, hLamBarFV3PReMassReCent_, "grLamBarFV3P");

    // Polar cosTheta^2 -- acceptance term
    GrData grLambdaP2Cos2(bpPb, hLambdaP2Cos2ReMassReCent_, "grLambdaP2Cos2");
    GrData grLambdaM2Cos2(bpPb, hLambdaM2Cos2ReMassReCent_, "grLambdaM2Cos2");
    GrData grLamBarP2Cos2(bpPb, hLamBarP2Cos2ReMassReCent_, "grLamBarP2Cos2");
    GrData grLamBarM2Cos2(bpPb, hLamBarM2Cos2ReMassReCent_, "grLamBarM2Cos2");
    GrData grLambdaF2Cos2(bpPb, hLambdaF2Cos2ReMassReCent_, "grLambdaF2Cos2");
    GrData grLamBarF2Cos2(bpPb, hLamBarF2Cos2ReMassReCent_, "grLamBarF2Cos2");

    // Polar P(M)-side wrt HFm(HFp)
    GrData grLambdaP2Cos    (bpPb, hLambdaP2CosReMassReCent_   , "grLambdaP2Cos");
    GrData grLambdaM2Cos    (bpPb, hLambdaM2CosReMassReCent_   , "grLambdaM2Cos");
    GrData grLamBarP2Cos    (bpPb, hLamBarP2CosReMassReCent_   , "grLamBarP2Cos");
    GrData grLamBarM2Cos    (bpPb, hLamBarM2CosReMassReCent_   , "grLamBarM2Cos");

    GrData grLambdaP3Cos    (bpPb, hLambdaP3CosReMassReCent_   , "grLambdaP3Cos");
    GrData grLambdaM3Cos    (bpPb, hLambdaM3CosReMassReCent_   , "grLambdaM3Cos");
    GrData grLamBarP3Cos    (bpPb, hLamBarP3CosReMassReCent_   , "grLamBarP3Cos");
    GrData grLamBarM3Cos    (bpPb, hLamBarM3CosReMassReCent_   , "grLamBarM3Cos");

    //// full eta
    GrData grLambdaF2MCos   (bpPb, hLambdaF2MCosReMassReCent_  , "grLambdaF2MCos");
    GrData grLambdaF3MCos   (bpPb, hLambdaF3MCosReMassReCent_  , "grLambdaF3MCos");
    GrData grLambdaF2PCos   (bpPb, hLambdaF2PCosReMassReCent_  , "grLambdaF2PCos");
    GrData grLambdaF3PCos   (bpPb, hLambdaF3PCosReMassReCent_  , "grLambdaF3PCos");

    GrData grLamBarF2MCos   (bpPb, hLamBarF2MCosReMassReCent_  , "grLamBarF2MCos");
    GrData grLamBarF3MCos   (bpPb, hLamBarF3MCosReMassReCent_  , "grLamBarF3MCos");
    GrData grLamBarF2PCos   (bpPb, hLamBarF2PCosReMassReCent_  , "grLamBarF2PCos");
    GrData grLamBarF3PCos   (bpPb, hLamBarF3PCosReMassReCent_  , "grLamBarF3PCos");

    // Polar SP
    GrData grLambdaP2CosSP  (bpPb, hLambdaP2CosSPReMassReCent_ , "grLambdaP2CosSP");
    GrData grLambdaM2CosSP  (bpPb, hLambdaM2CosSPReMassReCent_ , "grLambdaM2CosSP");
    GrData grLamBarP2CosSP  (bpPb, hLamBarP2CosSPReMassReCent_ , "grLamBarP2CosSP");
    GrData grLamBarM2CosSP  (bpPb, hLamBarM2CosSPReMassReCent_ , "grLamBarM2CosSP");

    GrData grLambdaP3CosSP  (bpPb, hLambdaP3CosSPReMassReCent_ , "grLambdaP3CosSP");
    GrData grLambdaM3CosSP  (bpPb, hLambdaM3CosSPReMassReCent_ , "grLambdaM3CosSP");
    GrData grLamBarP3CosSP  (bpPb, hLamBarP3CosSPReMassReCent_ , "grLamBarP3CosSP");
    GrData grLamBarM3CosSP  (bpPb, hLamBarM3CosSPReMassReCent_ , "grLamBarM3CosSP");

    //// full eta
    GrData grLambdaF2MCosSP (bpPb, hLambdaF2MCosSPReMassReCent_, "grLambdaF2MCosSP");
    GrData grLambdaF2PCosSP (bpPb, hLambdaF2PCosSPReMassReCent_, "grLambdaF2PCosSP");
    GrData grLamBarF2MCosSP (bpPb, hLamBarF2MCosSPReMassReCent_, "grLamBarF2MCosSP");
    GrData grLamBarF2PCosSP (bpPb, hLamBarF2PCosSPReMassReCent_, "grLamBarF2PCosSP");

    GrData grLambdaF3MCosSP (bpPb, hLambdaF3MCosSPReMassReCent_, "grLambdaF3MCosSP");
    GrData grLambdaF3PCosSP (bpPb, hLambdaF3PCosSPReMassReCent_, "grLambdaF3PCosSP");
    GrData grLamBarF3MCosSP (bpPb, hLamBarF3MCosSPReMassReCent_, "grLamBarF3MCosSP");
    GrData grLamBarF3PCosSP (bpPb, hLamBarF3PCosSPReMassReCent_, "grLamBarF3PCosSP");

    // Write Vn
    TFile * fsave = new TFile(Form("%s/graph.root", ftxt[s1]), "recreate");

    // Get Vn
    grLambdaPV2M.Write(fsave);
    grLambdaMV2P.Write(fsave);
    grLamBarPV2M.Write(fsave);
    grLamBarMV2P.Write(fsave);
    grLambdaFV2M.Write(fsave);
    grLambdaFV2P.Write(fsave);
    grLamBarFV2M.Write(fsave);
    grLamBarFV2P.Write(fsave);

    grLambdaPV3M.Write(fsave);
    grLambdaMV3P.Write(fsave);
    grLamBarPV3M.Write(fsave);
    grLamBarMV3P.Write(fsave);
    grLambdaFV3M.Write(fsave);
    grLambdaFV3P.Write(fsave);
    grLamBarFV3M.Write(fsave);
    grLamBarFV3P.Write(fsave);

    // Polar cosTheta^2 -- acceptance term
    grLambdaP2Cos2.Write(fsave);
    grLambdaM2Cos2.Write(fsave);
    grLamBarP2Cos2.Write(fsave);
    grLamBarM2Cos2.Write(fsave);
    grLambdaF2Cos2.Write(fsave);
    grLamBarF2Cos2.Write(fsave);

    // Polar P(M)-side wrt HFm(HFp)
    grLambdaP2Cos    .Write(fsave);
    grLambdaM2Cos    .Write(fsave);
    grLamBarP2Cos    .Write(fsave);
    grLamBarM2Cos    .Write(fsave);

    grLambdaP3Cos    .Write(fsave);
    grLambdaM3Cos    .Write(fsave);
    grLamBarP3Cos    .Write(fsave);
    grLamBarM3Cos    .Write(fsave);

    //// full eta
    grLambdaF2MCos   .Write(fsave);
    grLambdaF3MCos   .Write(fsave);
    grLambdaF2PCos   .Write(fsave);
    grLambdaF3PCos   .Write(fsave);

    grLamBarF2MCos   .Write(fsave);
    grLamBarF3MCos   .Write(fsave);
    grLamBarF2PCos   .Write(fsave);
    grLamBarF3PCos   .Write(fsave);

    // Polar SP
    grLambdaP2CosSP  .Write(fsave);
    grLambdaM2CosSP  .Write(fsave);
    grLamBarP2CosSP  .Write(fsave);
    grLamBarM2CosSP  .Write(fsave);

    grLambdaP3CosSP  .Write(fsave);
    grLambdaM3CosSP  .Write(fsave);
    grLamBarP3CosSP  .Write(fsave);
    grLamBarM3CosSP  .Write(fsave);

    //// full eta
    grLambdaF2MCosSP .Write(fsave);
    grLambdaF2PCosSP .Write(fsave);
    grLamBarF2MCosSP .Write(fsave);
    grLamBarF2PCosSP .Write(fsave);

    grLambdaF3MCosSP .Write(fsave);
    grLambdaF3PCosSP .Write(fsave);
    grLamBarF3MCosSP .Write(fsave);
    grLamBarF3PCosSP .Write(fsave);


    fsave->mkdir("LambdaMass");
    fsave->cd("/LambdaMass");
    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaMass_ [cent][ipt]->Write(Form("hLambdaMass_%i_%i",  cent, ipt));
            hLambdaMassP_[cent][ipt]->Write(Form("hLambdaMassP_%i_%i", cent, ipt));
            hLambdaMassM_[cent][ipt]->Write(Form("hLambdaMassM_%i_%i", cent, ipt));
        }
    }

    fsave->mkdir("LamBarMass");
    fsave->cd("/LamBarMass");
    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLamBarMass_ [cent][ipt]->Write(Form("hLamBarMass_%i_%i",  cent, ipt));
            hLamBarMassP_[cent][ipt]->Write(Form("hLamBarMassP_%i_%i", cent, ipt));
            hLamBarMassM_[cent][ipt]->Write(Form("hLamBarMassM_%i_%i", cent, ipt));
        }
    }

}

