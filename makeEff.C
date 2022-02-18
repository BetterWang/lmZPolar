
void makeEff(string s = "EPOS_pPb_V0_eff.root", string fout = "LmEPOS_pPb_eff.root")
{
    TFile * f = new TFile(s.c_str());

    TH3D * LamEff_num   = (TH3D*) f->Get("V0valid/Lm/LamEff_num");
    TH3D * LamEff_acc   = (TH3D*) f->Get("V0valid/Lm/LamEff_acc");
    TH3D * LamEff_denom = (TH3D*) f->Get("V0valid/Lm/LamEff_denom");
    TH3D * LamFak_num   = (TH3D*) f->Get("V0valid/Lm/LamFake_num");
    TH3D * LamFak_denom = (TH3D*) f->Get("V0valid/Lm/LamFake_denom");

    int mass_start = 64;
    int mass_end   = 82;

    const double * pt_bins   = LamEff_num->GetXaxis()->GetXbins()->GetArray();
    const double * eta_bins  = LamEff_num->GetYaxis()->GetXbins()->GetArray();
    const double * mass_bins = LamEff_num->GetZaxis()->GetXbins()->GetArray();

    int Npt_bins   = LamEff_num->GetXaxis()->GetNbins();
    int Neta_bins  = LamEff_num->GetYaxis()->GetNbins();
    int Nmass_bins = LamEff_num->GetZaxis()->GetNbins();

    auto pt_axis   = LamEff_num->GetXaxis();
    auto eta_axis  = LamEff_num->GetYaxis();
    auto mass_axis = LamEff_num->GetZaxis();

    double eff_eta_bins[11] = { -2.5, -2.0, -1.5, -1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5};
// pt                            0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12
//    double pTbin[14] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};

    TH1D* Eff_Num_eta  = new TH1D("Eff_Num_eta", "Eff_Num_eta", 10, eff_eta_bins); // 1.8 - 8.5
    TH1D* Eff_Acc_eta  = new TH1D("Eff_Acc_eta", "Eff_Acc_eta", 10, eff_eta_bins); // 1.8 - 8.5
    TH1D* Eff_Denom_eta= new TH1D("Eff_Denom_eta", "Eff_Denom_eta", 10, eff_eta_bins); // 1.8 - 8.5

    TH1D* Eff_Num_pT12  = new TH1D("Eff_Num_pT12", "Eff_Num_pT12", Npt_bins, pt_bins);
    TH1D* Eff_Acc_pT12  = new TH1D("Eff_Acc_pT12", "Eff_Acc_pT12", Npt_bins, pt_bins);
    TH1D* Eff_Denom_pT12= new TH1D("Eff_Denom_pT12", "Eff_Denom_pT12", Npt_bins, pt_bins);

    TH1D* Fak_Num_eta  = new TH1D("Fak_Num_eta", "Fak_Num_eta", 10, eff_eta_bins);
    TH1D* Fak_Denom_eta= new TH1D("Fak_Denom_eta", "Fak_Denom_eta", 10, eff_eta_bins);

    TH1D* Fak_Num_pT12  = new TH1D("Fak_Num_pT12", "Fak_Num_pT12", Npt_bins, pt_bins);
    TH1D* Fak_Denom_pT12= new TH1D("Fak_Denom_pT12", "Fak_Denom_pT12", Npt_bins, pt_bins);

    // eta range: -1.0 - 1.0, [16, 36)
    for ( int ieta = 15; ieta < 35; ieta++ ) {
        double eta = eta_axis->GetBinCenter(ieta+1);
        for ( int ipt = 0; ipt < Npt_bins; ipt++ ) {
            double pt = pt_axis->GetBinCenter(ipt+1);
            for ( int imass = mass_start; imass < mass_end; imass++ ) {
                double nLmEff = LamEff_num->GetBinContent( ipt+1, ieta+1, imass );
                double nLmAcc = LamEff_acc->GetBinContent( ipt+1, ieta+1, imass );
                double nLmDenom = LamEff_denom->GetBinContent( ipt+1, ieta+1, imass );

                Eff_Num_pT12->Fill( pt, nLmEff );
                Eff_Acc_pT12->Fill( pt, nLmAcc );
                Eff_Denom_pT12->Fill( pt, nLmDenom );

                double nLmFak = LamFak_num->GetBinContent( ipt+1, ieta+1, imass );
                double nLmFakDenom = LamFak_denom->GetBinContent( ipt+1, ieta+1, imass );

                Fak_Num_pT12->Fill( pt, nLmFak);
                Fak_Denom_pT12->Fill( pt, nLmFakDenom);
            }
        }
    }

    // pt range: 1.0 - 8.5  [4, 14)
    int pt_start = 4;
    int pt_end   = 13;

    for ( int ieta = 0; ieta < Neta_bins; ieta++ ) {
        double eta = eta_axis->GetBinCenter(ieta+1);
        for ( int ipt = pt_start; ipt < pt_end; ipt++ ) {
            double pt = pt_axis->GetBinCenter(ipt+1);
            for ( int imass = mass_start; imass < mass_end; imass++ ) {
                double nLmEff = LamEff_num->GetBinContent( ipt+1, ieta+1, imass );
                double nLmAcc = LamEff_acc->GetBinContent( ipt+1, ieta+1, imass );
                double nLmDenom = LamEff_denom->GetBinContent( ipt+1, ieta+1, imass );

                Eff_Num_eta->Fill( eta, nLmEff );
                Eff_Acc_eta->Fill( eta, nLmAcc );
                Eff_Denom_eta->Fill( eta, nLmDenom );

                double nLmFak = LamFak_num->GetBinContent( ipt+1, ieta+1, imass );
                double nLmFakDenom = LamFak_denom->GetBinContent( ipt+1, ieta+1, imass );

                Fak_Num_eta->Fill(eta, nLmFak);
                Fak_Denom_eta->Fill(eta, nLmFakDenom);
            }
        }
    }

    TFile * fsave = new TFile(fout.c_str(), "recreate");
    Eff_Num_eta  ->Write();
    Eff_Acc_eta  ->Write();
    Eff_Denom_eta->Write();

    Eff_Num_pT12  ->Write();
    Eff_Acc_pT12  ->Write();
    Eff_Denom_pT12->Write();

    Fak_Num_eta   ->Write();
    Fak_Denom_eta ->Write();

    Fak_Num_pT12  ->Write();
    Fak_Denom_pT12->Write();

    TH1D * hAcc_eta = (TH1D*) Eff_Acc_eta->Clone("hAcc_eta");
    TH1D * hEff_eta = (TH1D*) Eff_Num_eta->Clone("hEff_eta");
    TH1D * hREC_eta = (TH1D*) Eff_Num_eta->Clone("hREC_eta");

    hAcc_eta->Divide(Eff_Denom_eta);
    hEff_eta->Divide(Eff_Denom_eta);
    hREC_eta->Divide(  Eff_Acc_eta);

    TH1D * hAcc_pT12 = (TH1D*) Eff_Acc_pT12->Clone("hAcc_pT12");
    TH1D * hEff_pT12 = (TH1D*) Eff_Num_pT12->Clone("hEff_pT12");
    TH1D * hREC_pT12 = (TH1D*) Eff_Num_pT12->Clone("hREC_pT12");

    hAcc_pT12->Divide(Eff_Denom_pT12);
    hEff_pT12->Divide(Eff_Denom_pT12);
    hREC_pT12->Divide(  Eff_Acc_pT12);

    hAcc_eta->Write();
    hEff_eta->Write();
    hREC_eta->Write();

    hAcc_pT12->Write();
    hEff_pT12->Write();
    hREC_pT12->Write();

}
