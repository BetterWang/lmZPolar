#include "label.h"
#include "const.h"
#include "helper.h"

void plot(int s1)
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

    TFile * f = new TFile(Form("%s/graph.root", ftxt[s1]));
    // Get Vn
    GrData grLambdaPV2M(bpPb, f, "grLambdaPV2M");
    GrData grLambdaMV2P(bpPb, f, "grLambdaMV2P");
    GrData grLamBarPV2M(bpPb, f, "grLamBarPV2M");
    GrData grLamBarMV2P(bpPb, f, "grLamBarMV2P");
    GrData grLambdaFV2M(bpPb, f, "grLambdaFV2M");
    GrData grLambdaFV2P(bpPb, f, "grLambdaFV2P");
    GrData grLamBarFV2M(bpPb, f, "grLamBarFV2M");
    GrData grLamBarFV2P(bpPb, f, "grLamBarFV2P");

    GrData grLambdaPV3M(bpPb, f, "grLambdaPV3M");
    GrData grLambdaMV3P(bpPb, f, "grLambdaMV3P");
    GrData grLamBarPV3M(bpPb, f, "grLamBarPV3M");
    GrData grLamBarMV3P(bpPb, f, "grLamBarMV3P");
    GrData grLambdaFV3M(bpPb, f, "grLambdaFV3M");
    GrData grLambdaFV3P(bpPb, f, "grLambdaFV3P");
    GrData grLamBarFV3M(bpPb, f, "grLamBarFV3M");
    GrData grLamBarFV3P(bpPb, f, "grLamBarFV3P");

    // Polar cosTheta^2 -- acceptance term
    GrData grLambdaP2Cos2(bpPb, f, "grLambdaP2Cos2");
    GrData grLambdaM2Cos2(bpPb, f, "grLambdaM2Cos2");
    GrData grLamBarP2Cos2(bpPb, f, "grLamBarP2Cos2");
    GrData grLamBarM2Cos2(bpPb, f, "grLamBarM2Cos2");
    GrData grLambdaF2Cos2(bpPb, f, "grLambdaF2Cos2");
    GrData grLamBarF2Cos2(bpPb, f, "grLamBarF2Cos2");

    // Polar P(M)-side wrt HFm(HFp)
    GrData grLambdaP2Cos    (bpPb, f, "grLambdaP2Cos");
    GrData grLambdaM2Cos    (bpPb, f, "grLambdaM2Cos");
    GrData grLamBarP2Cos    (bpPb, f, "grLamBarP2Cos");
    GrData grLamBarM2Cos    (bpPb, f, "grLamBarM2Cos");

    GrData grLambdaP3Cos    (bpPb, f, "grLambdaP3Cos");
    GrData grLambdaM3Cos    (bpPb, f, "grLambdaM3Cos");
    GrData grLamBarP3Cos    (bpPb, f, "grLamBarP3Cos");
    GrData grLamBarM3Cos    (bpPb, f, "grLamBarM3Cos");

    //// full eta
    GrData grLambdaF2MCos   (bpPb, f, "grLambdaF2MCos");
    GrData grLambdaF3MCos   (bpPb, f, "grLambdaF3MCos");
    GrData grLambdaF2PCos   (bpPb, f, "grLambdaF2PCos");
    GrData grLambdaF3PCos   (bpPb, f, "grLambdaF3PCos");

    GrData grLamBarF2MCos   (bpPb, f, "grLamBarF2MCos");
    GrData grLamBarF3MCos   (bpPb, f, "grLamBarF3MCos");
    GrData grLamBarF2PCos   (bpPb, f, "grLamBarF2PCos");
    GrData grLamBarF3PCos   (bpPb, f, "grLamBarF3PCos");

    // Polar SP
    GrData grLambdaP2CosSP  (bpPb, f, "grLambdaP2CosSP");
    GrData grLambdaM2CosSP  (bpPb, f, "grLambdaM2CosSP");
    GrData grLamBarP2CosSP  (bpPb, f, "grLamBarP2CosSP");
    GrData grLamBarM2CosSP  (bpPb, f, "grLamBarM2CosSP");

    GrData grLambdaP3CosSP  (bpPb, f, "grLambdaP3CosSP");
    GrData grLambdaM3CosSP  (bpPb, f, "grLambdaM3CosSP");
    GrData grLamBarP3CosSP  (bpPb, f, "grLamBarP3CosSP");
    GrData grLamBarM3CosSP  (bpPb, f, "grLamBarM3CosSP");

    //// full eta
    GrData grLambdaF2MCosSP (bpPb, f, "grLambdaF2MCosSP");
    GrData grLambdaF2PCosSP (bpPb, f, "grLambdaF2PCosSP");
    GrData grLamBarF2MCosSP (bpPb, f, "grLamBarF2MCosSP");
    GrData grLamBarF2PCosSP (bpPb, f, "grLamBarF2PCosSP");

    GrData grLambdaF3MCosSP (bpPb, f, "grLambdaF3MCosSP");
    GrData grLambdaF3PCosSP (bpPb, f, "grLambdaF3PCosSP");
    GrData grLamBarF3MCosSP (bpPb, f, "grLamBarF3MCosSP");
    GrData grLamBarF3PCosSP (bpPb, f, "grLamBarF3PCosSP");

    // Plot
//    grLambdaPV2P.GrDump( ftxt[s1] );
//    grLambdaPV2M.GrDump( ftxt[s1] );
    grLambdaFV2P.GrDump( ftxt[s1], "vn" );
    grLambdaFV2M.GrDump( ftxt[s1], "vn" );
    grLambdaFV3P.GrDump( ftxt[s1], "vn" );
    grLambdaFV3M.GrDump( ftxt[s1], "vn" );

    // Polar cosTheta^2 -- acceptance term
    grLambdaP2Cos2.GrDump( ftxt[s1], "cos2" );
    grLambdaM2Cos2.GrDump( ftxt[s1], "cos2" );
    grLamBarP2Cos2.GrDump( ftxt[s1], "cos2" );
    grLamBarM2Cos2.GrDump( ftxt[s1], "cos2" );
    grLambdaF2Cos2.GrDump( ftxt[s1], "cos2" );
    grLamBarF2Cos2.GrDump( ftxt[s1], "cos2" );
}
