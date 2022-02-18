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

    TFile * f = new TFile(Form("%s/bGet_%i.root", ftxt[s1], 10));

    // Get
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

    // New
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

    for ( int cent = bpPb?6:0; cent < NCent; cent++ ) {
        for ( int ipt = 0; ipt < NpT; ipt++ ) {
            hLambdaPV2M_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaPV2M_%i_%i", cent, ipt));
            hLambdaMV2P_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaMV2P_%i_%i", cent, ipt));
            hLambdaFV2M_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaFV2M_%i_%i", cent, ipt));
            hLambdaFV2P_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaFV2P_%i_%i", cent, ipt));
            hLamBarPV2M_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarPV2M_%i_%i", cent, ipt));
            hLamBarMV2P_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarMV2P_%i_%i", cent, ipt));
            hLamBarFV2M_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarFV2M_%i_%i", cent, ipt));
            hLamBarFV2P_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarFV2P_%i_%i", cent, ipt));

            hLambdaPV3M_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaPV3M_%i_%i", cent, ipt));
            hLambdaMV3P_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaMV3P_%i_%i", cent, ipt));
            hLambdaFV3M_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaFV3M_%i_%i", cent, ipt));
            hLambdaFV3P_ [cent][ipt] = (TH1D*) f->Get(Form("hLambdaFV3P_%i_%i", cent, ipt));
            hLamBarPV3M_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarPV3M_%i_%i", cent, ipt));
            hLamBarMV3P_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarMV3P_%i_%i", cent, ipt));
            hLamBarFV3M_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarFV3M_%i_%i", cent, ipt));
            hLamBarFV3P_ [cent][ipt] = (TH1D*) f->Get(Form("hLamBarFV3P_%i_%i", cent, ipt));

            // ReMass
            hLambdaPV2MReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaPV2MReMass_%i_%i", cent, ipt), hLambdaPV2M_ [cent][ipt]);
            hLambdaMV2PReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaMV2PReMass_%i_%i", cent, ipt), hLambdaMV2P_ [cent][ipt]);
            hLambdaFV2MReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaFV2MReMass_%i_%i", cent, ipt), hLambdaFV2M_ [cent][ipt]);
            hLambdaFV2PReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaFV2PReMass_%i_%i", cent, ipt), hLambdaFV2P_ [cent][ipt]);
            hLamBarPV2MReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarPV2MReMass_%i_%i", cent, ipt), hLamBarPV2M_ [cent][ipt]);
            hLamBarMV2PReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarMV2PReMass_%i_%i", cent, ipt), hLamBarMV2P_ [cent][ipt]);
            hLamBarFV2MReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarFV2MReMass_%i_%i", cent, ipt), hLamBarFV2M_ [cent][ipt]);
            hLamBarFV2PReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarFV2PReMass_%i_%i", cent, ipt), hLamBarFV2P_ [cent][ipt]);

            hLambdaPV3MReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaPV3MReMass_%i_%i", cent, ipt), hLambdaPV3M_ [cent][ipt]);
            hLambdaMV3PReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaMV3PReMass_%i_%i", cent, ipt), hLambdaMV3P_ [cent][ipt]);
            hLambdaFV3MReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaFV3MReMass_%i_%i", cent, ipt), hLambdaFV3M_ [cent][ipt]);
            hLambdaFV3PReMass_ [cent][ipt] = RebinLmMass(Form("hLambdaFV3PReMass_%i_%i", cent, ipt), hLambdaFV3P_ [cent][ipt]);
            hLamBarPV3MReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarPV3MReMass_%i_%i", cent, ipt), hLamBarPV3M_ [cent][ipt]);
            hLamBarMV3PReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarMV3PReMass_%i_%i", cent, ipt), hLamBarMV3P_ [cent][ipt]);
            hLamBarFV3MReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarFV3MReMass_%i_%i", cent, ipt), hLamBarFV3M_ [cent][ipt]);
            hLamBarFV3PReMass_ [cent][ipt] = RebinLmMass(Form("hLamBarFV3PReMass_%i_%i", cent, ipt), hLamBarFV3P_ [cent][ipt]);
        }
    }

    // ReCent
    if ( bpPb ) {
        for ( int cent = 6; cent < NCent; cent++ ) {
            for ( int ipt = 0; ipt < NpT; ipt++ ) {
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

    // Write Vn
    TFile * fsave = new TFile(Form("%s/graph.root", ftxt[s1]), "recreate");

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

}

