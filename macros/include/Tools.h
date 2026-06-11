#pragma once

#include <cmath>
#include <iostream>
#include "TString.h"
#include "TFile.h"

class AnalysisTool {

public:

  // Constructors, destructor
  AnalysisTool(bool do_puppi_);
  AnalysisTool(const AnalysisTool &) = default;
  AnalysisTool & operator = (const AnalysisTool &) = default;
  ~AnalysisTool() = default;

  // Main functions
  void CalculateReconstructionQuality();
  void CalculateChi2Values();
  void CompareCHSPuppi();
  void PlotLimits(bool draw_data = false);
  void PlotPostfitDistributions();
  void PlotPostfitParameters();
  void ProduceThetaHistograms();

private:
  bool do_puppi;
  TString tag;
  TString path;
  std::vector<TString> CMcuts;
  std::vector<TString> CMcuts_str;

};
