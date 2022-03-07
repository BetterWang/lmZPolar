#include <TColor.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include "../../style.h"

#include "const.h"

TColor *red = new TColor(3001, 1, 0, 0, "red", 0.4);
TColor *green = new TColor(3002, 0, 0.5, 0, "green", 0.4);
TColor *blue = new TColor(3003, 0, 0, 1, "blue", 0.4);
TColor *gray = new TColor(3005, 0, 0, 0, "black", 0.4);


void setGr(TGraphErrors* gr, int kstyle, int kcolor, int ksz = 2)
{
    gr->SetMarkerStyle(kstyle);
    gr->SetMarkerColor(kcolor);
    gr->SetLineColor(kcolor);
    gr->SetMarkerSize(ksz);
}

void GrFilter(TGraphErrors* gr)
{
    for ( int i = 0; i < gr->GetN(); i++ ) {
        if ( TMath::IsNaN(gr->GetY()[i]) ) {
            gr->GetY()[i] = 0.;
            gr->GetEY()[i] = 0.;
        }
    }
    return;
}

struct GrData {
    bool          bpPb;
    TGraphErrors* grSig  [10] = {};
    TGraphErrors* grSBNeg[10] = {};
    TGraphErrors* grSBPos[10] = {};
    TGraphErrors* grSB   [10] = {};
    TString       Name;

    void Write(TFile * f) {
        f->cd();
        for ( int c = bpPb?6:0; c < (bpPb?NCentpPb2016:NCentPbPb2018Rebin); c++ ) {
//            GrFilter( grSig  [c] );
//            GrFilter( grSBNeg[c] );
//            GrFilter( grSBPos[c] );
//            GrFilter( grSB   [c] );
            grSig  [c]->Write(Form("%sSig_%i", Name.Data(), c));
            grSBNeg[c]->Write(Form("%sSBNeg_%i", Name.Data(), c));
            grSBPos[c]->Write(Form("%sSBPos_%i", Name.Data(), c));
            grSB   [c]->Write(Form("%sSB_%i", Name.Data(), c));
        }
    };

    void GrDump(TString dir = "", TString option = "vn")
    {
        TH2D* hframe = nullptr;
        if ( option == "vn" ) {
            hframe = new TH2D("hframe", "hframe", 1, 0, 10, 1, -0.1, 0.35);
            InitHist(hframe, "p_{T} (GeV/c)", "v_{n}");
        }
        if ( option == "cos2" ) {
            hframe = new TH2D("hframe", "hframe", 1, 0, 10, 1, -0.1, 0.55);
            InitHist(hframe, "p_{T} (GeV/c)", "<cos(#theta)^2>");
        }
        TCanvas * cT = MakeCanvas("cT", "cT", 600, 500);
        for ( int c = bpPb?6:0; c < (bpPb?NCentpPb2016:NCentPbPb2018Rebin); c++ ) {
            setGr( grSig[c], kFullSquare, kBlue );
            setGr( grSB [c], kOpenSquare, kRed  );
            hframe->Draw();
            grSig[c]->Draw("psame");
            grSB [c]->Draw("psame");
            cT->SaveAs( (dir + "/" + Name + c + ".pdf" ).Data() );
        }
        delete hframe;
        delete cT;
    };

    GrData(bool b, TFile *f, TString name) :
        bpPb(b),
        Name(name)
    {
        f->cd();
        for ( int c = (bpPb?6:0); c < (bpPb?NCentpPb2016:NCentPbPb2018Rebin); c++ ) {
            grSig  [c] = (TGraphErrors*) f->Get(Form("%sSig_%i", Name.Data(), c));
            grSBNeg[c] = (TGraphErrors*) f->Get(Form("%sSBNeg_%i", Name.Data(), c));
            grSBPos[c] = (TGraphErrors*) f->Get(Form("%sSBPos_%i", Name.Data(), c));
            grSB   [c] = (TGraphErrors*) f->Get(Form("%sSB_%i", Name.Data(), c));
        }
    };

    GrData(bool b, vector<vector<TH1D*>> &h, TString name) :
        bpPb(b),
        Name(name)
    {
        for ( int c = (bpPb?6:0); c < (bpPb?NCentpPb2016:NCentPbPb2018Rebin); c++ ) {
            grSig  [c] = new TGraphErrors(NpT);
            grSBNeg[c] = new TGraphErrors(NpT);
            grSBPos[c] = new TGraphErrors(NpT);
            grSB   [c] = new TGraphErrors(NpT);
            for ( int ipt = 0; ipt < NpT; ipt++ ) {
                if ( bpPb ) {
                    grSig  [c]->GetX()[ipt] = pTX_pPb[c][ipt];
                    grSBNeg[c]->GetX()[ipt] = pTX_pPb[c][ipt];
                    grSBPos[c]->GetX()[ipt] = pTX_pPb[c][ipt];
                    grSB   [c]->GetX()[ipt] = pTX_pPb[c][ipt];
                } else {
                    grSig  [c]->GetX()[ipt] = pTX_PbPb[c][ipt];
                    grSBNeg[c]->GetX()[ipt] = pTX_PbPb[c][ipt];
                    grSBPos[c]->GetX()[ipt] = pTX_PbPb[c][ipt];
                    grSB   [c]->GetX()[ipt] = pTX_PbPb[c][ipt];
                }

                grSBNeg[c]->GetY()[ipt] = h[c][ipt]->GetBinContent(1);
                grSig  [c]->GetY()[ipt] = h[c][ipt]->GetBinContent(3);
                grSBPos[c]->GetY()[ipt] = h[c][ipt]->GetBinContent(5);

                grSBNeg[c]->GetEY()[ipt] = h[c][ipt]->GetBinError(1);
                grSig  [c]->GetEY()[ipt] = h[c][ipt]->GetBinError(3);
                grSBPos[c]->GetEY()[ipt] = h[c][ipt]->GetBinError(5);
            }
        }
        for ( int c = bpPb?6:0; c < (bpPb?NCentpPb2016:NCentPbPb2018Rebin); c++ ) {
            for ( int i = 0; i < NpT; i++ ) {
                double pos_w2 = grSBPos[c]->GetEY()[i] * grSBPos[c]->GetEY()[i];
                double neg_w2 = grSBNeg[c]->GetEY()[i] * grSBNeg[c]->GetEY()[i];
                double sum = (grSBPos[c]->GetY()[i] / pos_w2 + grSBNeg[c]->GetY()[i] / neg_w2) / (1./pos_w2 + 1./neg_w2);
                if ( (pos_w2==0.) or (neg_w2==0.) ) {
                    grSB[c]->GetY()[i] = -999.;
                    grSB[c]->GetEY()[i] = 9999.;
                } else {
                    grSB[c]->GetY()[i] = sum;
                    grSB[c]->GetEY()[i]= 1. / sqrt( 1./pos_w2 + 1./neg_w2 );
                }
            }
        }
    };
};




// rebin mass histo
TH1D* RebinLmMass(const char* name, TH1D* h, TString options = "")
{
    double sum = 0.;
    double sum_w = 0.;
    TH1D * ret = new TH1D(name, "", 5, Lm_mass_ranges);
    // SBNeg
    for ( int i = 1; i <= Lm_SBNeg_Bin; i++ ) {
        double sigma2 = h->GetBinError(i) * h->GetBinError(i);
        sum += h->GetBinContent(i) / sigma2;
        sum_w += 1./sigma2;
    }
    if ( sum_w == 0. ) {
        ret->SetBinContent( 1, -999. );
        ret->SetBinError  ( 1, 9999. );
    } else {
        ret->SetBinContent( 1, sum / sum_w );
        ret->SetBinError  ( 1, 1./sqrt(sum_w) );
    }

    //MixNeg
    sum = 0.;
    sum_w = 0.;
    for ( int i = Lm_SBNeg_Bin+1; i < Lm_sig_min_Bin; i++ ) {
        double sigma2 = h->GetBinError(i) * h->GetBinError(i);
        sum += h->GetBinContent(i) / sigma2;
        sum_w += 1./sigma2;
    }
    if ( sum_w == 0. ) {
        ret->SetBinContent( 2, -999. );
        ret->SetBinError  ( 2, 9999. );
    } else {
        ret->SetBinContent( 2, sum / sum_w );
        ret->SetBinError  ( 2, 1./sqrt(sum_w) );
    }

    //Sig
    sum = 0.;
    sum_w = 0.;
    for ( int i = Lm_sig_min_Bin; i <= Lm_sig_max_Bin; i++ ) {
        double sigma2 = h->GetBinError(i) * h->GetBinError(i);
        sum += h->GetBinContent(i) / sigma2;
        sum_w += 1./sigma2;
    }
    if ( sum_w == 0. ) {
        ret->SetBinContent( 3, -999. );
        ret->SetBinError  ( 3, 9999. );
    } else {
        ret->SetBinContent( 3, sum / sum_w );
        ret->SetBinError  ( 3, 1./sqrt(sum_w) );
    }

    //MixPos
    sum = 0.;
    sum_w = 0.;
    for ( int i = Lm_sig_max_Bin+1; i < Lm_SBPos_Bin; i++ ) {
        double sigma2 = h->GetBinError(i) * h->GetBinError(i);
        sum += h->GetBinContent(i) / sigma2;
        sum_w += 1./sigma2;
    }
    if ( sum_w == 0. ) {
        ret->SetBinContent( 4, -999. );
        ret->SetBinError  ( 4, 9999. );
    } else {
        ret->SetBinContent( 4, sum / sum_w );
        ret->SetBinError  ( 4, 1./sqrt(sum_w) );
    }

    //SBPos
    sum = 0.;
    sum_w = 0.;
    for ( int i = Lm_SBPos_Bin; i <= 160; i++ ) {
        double sigma2 = h->GetBinError(i) * h->GetBinError(i);
        sum += h->GetBinContent(i) / sigma2;
        sum_w += 1./sigma2;
    }
    if ( sum_w == 0. ) {
        ret->SetBinContent( 5, -999. );
        ret->SetBinError  ( 5, 9999. );
    } else {
        ret->SetBinContent( 5, sum / sum_w );
        ret->SetBinError  ( 5, 1./sqrt(sum_w) );
    }

    return ret;
}


struct SaveToFile
{
    TFile * fsave;
    SaveToFile(TFile *f) { fsave = f; };
    void WriteGr( TGraphErrors * gr, string name ) {
        fsave->cd();
        gr->Write(name.c_str());
    };
    ~SaveToFile() {fsave->Close();};
};


void DropPoints( TGraphErrors* gr, int N )
{
    if ( N > 0 ) {
        for ( int i = 0; i < N; i++ ) {
            gr->RemovePoint(0);
        }
    } else {
        for ( int i = 0; i < -N; i++ ) {
            gr->RemovePoint(gr->GetN()-1);
        }
    }
}

TGraphErrors * getGr(TFile *f, string s, int kstyle, int kcolor, int ksz = 2)
{
    TGraphErrors * ret = (TGraphErrors*) f->Get(s.c_str());
    ret->SetMarkerStyle(kstyle);
    ret->SetMarkerColor(kcolor);
    ret->SetLineColor(kcolor);
    ret->SetMarkerSize(ksz);
    return ret;
}

TGraphErrors* getRatio(TGraphErrors* gr1, TGraphErrors* gr2, int option = 0)
    /*
     * 0 - errors are independent
     * 1 - use only gr1 error
     * 2 - use only gr2 error
     * 3 - errors are correlated
     */
{
    if ( gr1->GetN() != gr2->GetN() ) {
        cout << " --> Warning mismatch getRatio" << endl;
        return nullptr;
    }
    TGraphErrors * ret = new TGraphErrors(gr2->GetN());
    for ( int i = 0; i < gr2->GetN(); i++ ) {
        if ( gr1->GetX()[i] != gr2->GetX()[i] ) cout << "--> Warning mismatch X" << endl;
        ret->GetX()[i] = gr2->GetX()[i];
        ret->GetEX()[i] = gr2->GetEX()[i];
        ret->GetY()[i] = gr1->GetY()[i] / gr2->GetY()[i];
        double y1 = gr1->GetY()[i];
        double y2 = gr2->GetY()[i];
        double e1 = gr1->GetEY()[i];
        double e2 = gr2->GetEY()[i];
        if ( option == 1 ) {
            e2 = 0;
        } else if ( option == 2 ) {
            e1 = 0;
        }
        ret->GetEY()[i] = y1/y2*sqrt( e1*e1/y1/y1 + e2*e2/y2/y2 );
        if ( option == 3 ) {
            ret->GetEY()[i] = y1/y2*sqrt( abs(e1*e1/y1/y1 + e2*e2/y2/y2 - 2*e1*e1/y1/y2) );
        }
    }
    return ret;
}

