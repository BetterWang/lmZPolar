#include "label.h"

void verify(int s1) {
    int s3 = 10;
    cout << ftxt[s1] << endl;
    TFile *f = new TFile(Form("%s/output_%i.root", ftxt[s1], s3));
    if (!f->IsZombie()) {
        TH1D * h = (TH1D*) f->Get("output_10_10/hCent");
        unsigned int Nevt = h->GetEntries();
        cout << " s1 = " << s1 << "\tNevt = " << Nevt << endl;
    } else {
        cout << " s1 = " << s1 << "\tNot exist!!!" << endl;
    }
}
