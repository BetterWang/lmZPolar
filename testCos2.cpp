#include <TH1D.h>

void testCos2()
{
    double sum = 0;
    for ( int i = 0; i < 1000000; i++ ) {
        double theta = gRandom->Rndm() * TMath::Pi();
        sum += TMath::Cos(theta) * TMath::Cos(theta);
        cout << " sum = " << sum << "\t" << "avg = " << sum / (i+1) << endl;
    }
}
