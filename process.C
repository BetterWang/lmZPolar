#include "label.h"
#include "polar.h"
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

using namespace std;

void process(int s1 = 0, int s3 = 10, int sVz = 0)
{
    cout << " s1 = " << s1 << " s3 = " << s3 << endl;
    TH1::SetDefaultSumw2();
    addchain(s1);

    PolarEvent evt;

    chV->SetBranchAddress("EPOrg",      &evt.EPOrg);
    chV->SetBranchAddress("EPOrgSin",   &evt.EPOrgSin);
    chV->SetBranchAddress("EPOrgCos",   &evt.EPOrgCos);
    chV->SetBranchAddress("EPFlat",     &evt.EPFlat);
    chV->SetBranchAddress("EPFlatSin",  &evt.EPFlatSin);
    chV->SetBranchAddress("EPFlatCos",  &evt.EPFlatCos);
    chV->SetBranchAddress("pdgId",      &evt.pdgId);
    chV->SetBranchAddress("pTrkPt",     &evt.pTrkPt);
    chV->SetBranchAddress("pTrkPx",     &evt.pTrkPx);
    chV->SetBranchAddress("pTrkPy",     &evt.pTrkPy);
    chV->SetBranchAddress("pTrkPz",     &evt.pTrkPz);
    chV->SetBranchAddress("pTrkEta",    &evt.pTrkEta);
    chV->SetBranchAddress("nTrkPt",     &evt.nTrkPt);
    chV->SetBranchAddress("nTrkPx",     &evt.nTrkPx);
    chV->SetBranchAddress("nTrkPy",     &evt.nTrkPy);
    chV->SetBranchAddress("nTrkPz",     &evt.nTrkPz);
    chV->SetBranchAddress("nTrkEta",    &evt.nTrkEta);
    chV->SetBranchAddress("pPxCM",      &evt.pPxCM);
    chV->SetBranchAddress("pPyCM",      &evt.pPyCM);
    chV->SetBranchAddress("pPzCM",      &evt.pPzCM);
    chV->SetBranchAddress("nPxCM",      &evt.nPxCM);
    chV->SetBranchAddress("nPyCM",      &evt.nPyCM);
    chV->SetBranchAddress("nPzCM",      &evt.nPzCM);
    chV->SetBranchAddress("pt",         &evt.pt);
    chV->SetBranchAddress("phi",        &evt.phi);
    chV->SetBranchAddress("eta",        &evt.eta);
    chV->SetBranchAddress("rapidity",   &evt.rapidity);
    chV->SetBranchAddress("mass",       &evt.mass);
    chV->SetBranchAddress("Cent",       &evt.Cent);
    chV->SetBranchAddress("vz",         &evt.vz);
    chV->SetBranchAddress("RunId",      &evt.RunId);
    chV->SetBranchAddress("EventId",    &evt.EventId);
    chV->SetBranchAddress("Lumi",       &evt.Lumi);

    PolarData* data[s3+1];
    for ( int i = 0; i <= s3; i++ ) {
        data[i] = new PolarData(string("pd")+to_string(i));
    }

    int ievt = 0;
    while ( chV->GetEntry(ievt++) )
    {
        if ( !(ievt%1000000) ) cout << "!! ievt = " << ievt << endl;

        int s2 = ievt%s3;
        data[s2]->Fill(evt);
        int ret = data[s3]->Fill(evt);
        if ( ret == 1 ) {
            std::cout << " ---> Failed event in " << chV->GetFile()->GetName() << " RunId = " << evt.RunId << " Lumi = " << evt.Lumi << " EventId = " << evt.EventId << std::endl;
        }
    }

    TFile * fwrite = new TFile(Form("%s/output_%i.root", ftxt[s1], s3), "recreate");
    for ( int i = 0; i <= s3; i++ ) {
        auto dir = fwrite->mkdir(Form("output_%i_%i", i, s3));
        dir->cd();
        data[i]->Write();
    }
}
