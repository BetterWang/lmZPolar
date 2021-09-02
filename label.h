
char const * fname[] = {
    "test",         // 0 test
    "",             // 1 HIMB0
    "",             // 2 HIMB1
    "",             // 3 HIMB2
    "",             // 4 HIMB3
    "../PbPb2018/HIMinimumBias4/crab_HIMB4_LmPolarTree_v1/210529_031055/0000/",             // 5 HIMB4
};

char const * ftxt[] = {
    "txt/test/",    // 0
    "txt/HIMB0/",   // 1
    "txt/HIMB1/",   // 2
    "txt/HIMB2/",   // 3
    "txt/HIMB3/",   // 4
    "txt/HIMB4/",   // 5
};


#include <TChain.h>
#include <iostream>

TChain * chV = new TChain("trV");
void addchain(int s1)
{
    std::cout << fname[s1] << std::endl;
    chV->Add(Form("%s/*.root/tree/trV", fname[s1]));
}
