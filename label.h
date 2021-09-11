
char const * fname[] = {
    "test",         // 0 test
    "",             // 1 HIMB0
    "",             // 2 HIMB1
    "../PbPb2018/HIMinimumBias2/crab_HIMB2_LmPolarTree_v1/210911_162435/0000/",             // 3 HIMB2
    "",             // 4 HIMB3
    "../PbPb2018/HIMinimumBias4/crab_HIMB4_LmPolarTree_v1/210529_031055/0000/",             // 5 HIMB4
    "",             // 6 HIMB5
    "",             // 7 HIMB6
    "",             // 8 HIMB7
    "",             // 9 HIMB8
    "",             // 10 HIMB9
    "",             // 11 HIMB10
    "",             // 12 HIMB11
    "",             // 13 HIMB12
    "",             // 14 HIMB13
    "",             // 15 HIMB14
    "",             // 16 HIMB15
    "",             // 17 HIMB16
    "",             // 18 HIMB17
    "",             // 19 HIMB18
    "",             // 20 HIMB19
    "NA",           // 21 HIMB merge
    "../PbPb2018/HIMinimumBias2/crab_HIMB2_V0Tree_WrongSign_v9/210911_053022/0000/",           // 22 HIMB2 WrongSign
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
};


#include <TChain.h>
#include <iostream>

TChain * chV = new TChain("trV");
void addchain(int s1)
{
    std::cout << fname[s1] << std::endl;
    chV->Add(Form("%s/*.root/tree/trV", fname[s1]));
}
