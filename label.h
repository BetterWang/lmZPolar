
char const * fname[] = {
    "test",         // 0 test
    "",             // 1 HIMB0
    "../PbPb2018/HIMinimumBias1/crab_HIMB1_LmPolarTree_v1/210916_101315/0000/",             // 2 HIMB1
    "../PbPb2018/HIMinimumBias2/crab_HIMB2_LmPolarTree_v1/210911_162435/0000/",             // 3 HIMB2
    "../PbPb2018/HIMinimumBias3/crab_HIMB3_LmPolarTree_v1/210916_081822/0000/",             // 4 HIMB3
    "../PbPb2018/HIMinimumBias4/crab_HIMB4_LmPolarTree_v2/211007_091317/0000/",             // 5 HIMB4
    "../PbPb2018/HIMinimumBias5/crab_HIMB5_LmPolarTree_v1/210918_051058/0000/",             // 6 HIMB5
    "../PbPb2018/HIMinimumBias6/crab_HIMB6_LmPolarTree_v1/210927_133247/0000/",             // 7 HIMB6
    "../PbPb2018/HIMinimumBias7/crab_HIMB7_LmPolarTree_v1/210927_133414/0000/",             // 8 HIMB7
    "../PbPb2018/HIMinimumBias8/crab_HIMB8_LmPolarTree_v1/210927_133535/0000/",             // 9 HIMB8
    "../PbPb2018/HIMinimumBias9/crab_HIMB9_LmPolarTree_v1/210930_092403/0000/",             // 10 HIMB9
    "../PbPb2018/HIMinimumBias10/crab_HIMB10_LmPolarTree_v1/210930_112327/0000/",           // 11 HIMB10
    "../PbPb2018/HIMinimumBias11/crab_HIMB11_LmPolarTree_v1/210930_112509/0000/",             // 12 HIMB11
    "",             // 13 HIMB12
    "",             // 14 HIMB13
    "",             // 15 HIMB14
    "",             // 16 HIMB15
    "",             // 17 HIMB16
    "",             // 18 HIMB17
    "",             // 19 HIMB18
    "",             // 20 HIMB19
    "NA",           // 21 HIMB merge
    "../PbPb2018/HIMinimumBias2/crab_HIMB2_V0Tree_WrongSign_v2/211011_074510/0000/",           // 22 HIMB2 WrongSign
    "",           // 23
    "",           // 24
    "",           // 25
    "",           // 26
    "",           // 27
    "",           // 28
    "",           // 29
    "",           // 30 PA HM0 120
    "",           // 31 PA HM0 150
    "../PA_Polar/PAHighMultiplicity1/crab_PAHM1_LmPolarTree_HM185_v1/211115_221809/0000/",           // 32 PA HM1 185
    "../PA_Polar/PAHighMultiplicity2/crab_PAHM2_LmPolarTree_HM185_v1/211118_152421/0000/",           // 33 PA HM2 185
    "../PA_Polar/PAHighMultiplicity3/crab_PAHM3_LmPolarTree_HM185_v1/211118_190407/0000/",           // 34 PA HM3 185
    "../PA_Polar/PAHighMultiplicity4/crab_PAHM4_LmPolarTree_HM185_v1/211118_191658/0000/",           // 35 PA HM4 185
    "../PA_Polar/PAHighMultiplicity5/crab_PAHM5_LmPolarTree_HM185_v1/211118_193343/0000/",           // 36 PA HM5 185
    "../PA_Polar/PAHighMultiplicity6/crab_PAHM6_LmPolarTree_HM185_v1/211118_193559/0000/",           // 37 PA HM6 185
    "",           // 38 PA HM7 250
    "",           // 39 PA HM0 120 reverse
    "",           // 40 PA HM0 150 reverse
    "../PA_Polar/PAHighMultiplicity1/crab_PAHM1_LmPolarTree_HM185_reverse_v1/211115_224239/0000/",           // 41 PA HM1 185 reverse
    "../PA_Polar/PAHighMultiplicity2/crab_PAHM2_LmPolarTree_HM185_reverse_v1/211207_193924/0000/",           // 42 PA HM2 185 reverse
    "../PA_Polar/PAHighMultiplicity3/crab_PAHM3_LmPolarTree_HM185_reverse_v1/211208_201019/0000/",           // 43 PA HM3 185 reverse
    "",           // 44 PA HM4 185 reverse
    "",           // 45 PA HM5 185 reverse
    "",           // 46 PA HM6 185 reverse
    "",           // 47 PA HM7 250 reverse
};

char const * ftxt[] = {
    "txt/test/",    // 0
    "txt/HIMB0/",   // 1
    "txt/HIMB1/",   // 2
    "txt/HIMB2/",   // 3
    "txt/HIMB3/",   // 4
    "txt/HIMB4/",   // 5
    "txt/HIMB5/",   // 6
    "txt/HIMB6/",   // 7
    "txt/HIMB7/",   // 8
    "txt/HIMB8/",   // 9
    "txt/HIMB9/",   // 10
    "txt/HIMB10/",  // 11
    "txt/HIMB11/",  // 12
    "txt/HIMB12/",  // 13
    "txt/HIMB13/",  // 14
    "txt/HIMB14/",  // 15
    "txt/HIMB15/",  // 16
    "txt/HIMB16/",  // 17
    "txt/HIMB17/",  // 18
    "txt/HIMB18/",  // 19
    "txt/HIMB19/",  // 20
    "txt/HIMB_merge/",  // 21
    "txt/HIMB2_WS/",    // 22
    "",    // 23
    "",    // 24
    "",    // 25
    "",    // 26
    "",    // 27
    "",    // 28
    "",    // 29
    "txt/PAHM0_120",    // 30
    "txt/PAHM0_150",    // 31
    "txt/PAHM1_185",    // 32
    "txt/PAHM2_185",    // 33
    "txt/PAHM3_185",    // 34
    "txt/PAHM4_185",    // 35
    "txt/PAHM5_185",    // 36
    "txt/PAHM6_185",    // 37
    "txt/PAHM7_250",    // 38
    "txt/PAHM0_120_reverse",    // 39
    "txt/PAHM0_150_reverse",    // 40
    "txt/PAHM1_185_reverse",    // 41
    "txt/PAHM2_185_reverse",    // 42
    "txt/PAHM3_185_reverse",    // 43
    "txt/PAHM4_185_reverse",    // 44
    "txt/PAHM5_185_reverse",    // 45
    "txt/PAHM6_185_reverse",    // 46
    "txt/PAHM7_250_reverse",    // 47
};


#include <TChain.h>
#include <iostream>

TChain * chV = new TChain("trV");

bool bPbPb = false;
bool bpPb = false;
bool bpPbReverse = false;

void addchain(int s1)
{
    std::cout << fname[s1] << std::endl;
    chV->Add(Form("%s/*.root/tree/trV", fname[s1]));
    TString tt(fname[s1]);
    if ( tt.Contains("PbPb2018") ) bPbPb = true;
    if ( tt.Contains("PA_Polar") ) {
        bpPb = true;
        if ( tt.Contains("reverse") ) {
            bpPbReverse = true;
        } else {
            bpPbReverse = false;
        }
    }
}
