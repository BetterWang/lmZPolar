#ifndef __LmZPolarConst_
#define __LmZPolarConst_

// PbPb 2018
//                                     0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,   13,   14,   15,   16
//ptBin = cms.untracked.vdouble(0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5)
  const double pTX[] =              {0.3, 0.5, 0.7, 0.9, 1.2, 1.6, 2.0, 2.5, 3.2, 4.1, 5.3, 6.5, 7.7};
  const double ptBin[]=             {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};

  const int NpT = 13;

// PbPb 2018
// CentBin                      0,  1,   2,   3,   4,   5,   6,   7,   8,   9,   10,   11,   12,   13,   14,   15,   16,  17
  const Int_t CentPbPb2018[] = {0, 10,  20,  30,  40,  50,  60,  70,  80,  90,  100,  110,  120,  130,  140,  150,  160};
  const double CentPbPbX[]   = {  2.5, 7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5};

  const int NCentPbPb2018 = 16;

// mass window
  const double Lm_mass_min = 1.0800;
  const double Lm_mass_max = 1.1600;
  const double Lm_sig_min  = 1.1115;
  const double Lm_sig_max  = 1.1200;

#endif
