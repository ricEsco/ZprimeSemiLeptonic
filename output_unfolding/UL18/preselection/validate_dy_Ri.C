// validate_dy_Ri.C  — use: root -l -b -q 'validate_dy_Ri.C("uhh2.AnalysisModule.root")'
#include <TFile.h>
#include <TH1.h>
#include <cstdio>
#include <cmath>
void validate_dy_Ri(const char* file="uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL18.root"){
  const char* tags[6] = {"DY_0_350_General","DY_350_500_General","DY_500_750_General",
                         "DY_750_1000_General","DY_1000_1500_General","DY_1500Inf_General"};
  TFile f(file);
  printf("%-20s  %12s  %12s  %12s  %s\n","Bin","Nneg","Npos","R=N-/N+","SumOK");
  for (int i=0;i<6;++i){
    TH1* hE = (TH1*)f.Get((std::string(tags[i])+"/Events").c_str());
    TH1* hS = (TH1*)f.Get((std::string(tags[i])+"/dYgen").c_str());
    if(!hE || !hS){ printf("%-20s  MISSING\n", tags[i]); continue; }
    double Ne   = hE->GetBinContent(1);
    double Nneg = hS->GetBinContent(1);   // bin1 = Δ|y|<0
    double Npos = hS->GetBinContent(2);   // bin2 = Δ|y|≥0 (dY=0 goes here)
    double R    = (Npos>0) ? Nneg/Npos : 0.0;
    bool sumOK  = (std::fabs((Nneg+Npos)-Ne) <= 1e-6*std::max(1.0, std::fabs(Ne)));
    printf("%-20s  %12.6g  %12.6g  %12.8g  %s\n", tags[i], Nneg, Npos, R, sumOK?"OK":"check");
  }

  // Combine mapping template per bin i (Ac in percent; use (1±@1) if Ac is fractional)
  puts("\n# --PO map template per mass bin i (uses R_i = Nneg/Npos)");
  for (int i=0;i<6;++i){
    printf("--PO map='.*/TTbar_neg_M%d:r_neg_%d[1,0,50]' \\\n", i,i);
    printf("--PO map='.*/TTbar_pos_M%d:r_pos_%d=expr;;r_pos_%d(\"R_%d*@0*(100+@1)/(100-@1)\", r_neg_%d, Ac_%d[0,-5,5])' \\\n", i,i,i,i,i,i);
  }
  puts("# Replace R_i with the numbers printed above for each bin.");
}
