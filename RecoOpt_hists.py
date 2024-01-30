from ROOT import *
import numpy as np

# Choose whether to plot with log scale or not
# log = True
log = False

fileDir =           "/nfs/dust/cms/user/ricardo/analysis_output/RecoOpt/muon/workdir3_Analysis_UL18_muon_RecOpt_signal"
if log:
    plotDirectory = "/nfs/dust/cms/user/ricardo/analysis_output/RecoOpt/plots/UL18/muon/log"
else:
    plotDirectory = "/nfs/dust/cms/user/ricardo/analysis_output/RecoOpt/plots/UL18/muon"
txtfile =           "/nfs/dust/cms/user/ricardo/analysis_output/RecoOpt/EfficiencyInfo.txt"

print "---"
print "Input files will come from", fileDir
print "Output files will go into", plotDirectory
print "Integral and efficiency info will be saved to", txtfile
print "---"

### This array has filenames of samples that will contribute to each plot
sample_names = ["TTToSemiLeptonic", "TTToSemiLeptonic_2"]

### This array contains the names of the histograms to be plotted along with their : ["titles", number of bins, [range]]
hists = {   # "res_jet_bscore"                              : ["AK4 b-score (resolved)",                                          50, [ 0,   1]],
            # "mer_subjet_bscore"                           : ["AK8 subjet b-score (merged)",                                     50, [ 0,   1]],
            # "bscore_max"                                  : ["max b-scores (both)",                                             50, [ 0,   1]],
            # "res_jet_bscore_matchable"                    : ["AK4 b-score for matchable events (resolved)",                     50, [ 0,   1]],
            # "mer_subjet_bscore_matchable"                 : ["AK8 subjet b-score for matchable events (merged)",                50, [ 0,   1]],
            # "bscore_max_matchable"                        : ["max b-scores for matchable events (both)",                        50, [ 0,   1]],
    
            # "BestChi2"                                    : ["#chi^{2}",                                                        50, [ 0, 100]],
            # "BestChi2_resolved"                           : ["#chi^{2} (resolved)",                                             50, [ 0, 100]],
            # "BestChi2_merged"                             : ["#chi^{2} (merged)",                                               50, [ 0, 100]],
    
            # "Njets_chi2"                                  : ["N_{jets} of Best #Chi^{2} Candidate",                             10, [ 0,  10]],
            # "Njets_chi2_resolved"                         : ["N_{jets} of Best #Chi^{2} Candidate (resolved)",                  10, [ 0,  10]],
            # "Njets_chi2_merged"                           : ["N_{jets} of Best #Chi^{2} Candidate (merged)",                    10, [ 0,  10]],

            # "BestChi2_passedMatching"                     : ["#chi^{2} passed matching",                                        50, [ 0, 100]],
            # "BestChi2_passedMatching_resolved"            : ["#chi^{2} passed matching (resolved)",                             50, [ 0, 100]],
            # "BestChi2_passedMatching_merged"              : ["#chi^{2} passed matching (merged)",                               50, [ 0, 100]],

            # "Njets_chi2_passedChi2cut"                    : ["N_{jets} of Best #Chi^{2} Candidate passed #chi^{2} cut",             10, [ 0,  10]],
            # "Njets_chi2_passedChi2cut_resolved"           : ["N_{jets} of Best #Chi^{2} Candidate passed #chi^{2} cut (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_passedChi2cut_merged"             : ["N_{jets} of Best #Chi^{2} Candidate passed #chi^{2} cut (merged)",    10, [ 0,  10]],

            # "Njets_chi2_had_passedChi2cut"                : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg passed #chi^{2} cut",             10, [ 0,  10]],
            # "Njets_chi2_had_passedChi2cut_resolved"       : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg passed #chi^{2} cut (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_had_passedChi2cut_merged"         : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg passed #chi^{2} cut (merged)",    10, [ 0,  10]],

            # "Njets_chi2_lep_passedChi2cut"                : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg passed #chi^{2} cut",             10, [ 0,  10]],
            # "Njets_chi2_lep_passedChi2cut_resolved"       : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg passed #chi^{2} cut (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_lep_passedChi2cut_merged"         : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg passed #chi^{2} cut (merged)",    10, [ 0,  10]],

            # "Njets_chi2_failedChi2cut"                    : ["N_{jets} of Best #Chi^{2} Candidate failed #chi^{2} cut",             10, [ 0,  10]],
            # "Njets_chi2_failedChi2cut_resolved"           : ["N_{jets} of Best #Chi^{2} Candidate failed #chi^{2} cut (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_failedChi2cut_merged"             : ["N_{jets} of Best #Chi^{2} Candidate failed #chi^{2} cut (merged)",    10, [ 0,  10]],

            # "Njets_chi2_had_failedChi2cut"                : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg failed #chi^{2} cut",             10, [ 0,  10]],
            # "Njets_chi2_had_failedChi2cut_resolved"       : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg failed #chi^{2} cut (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_had_failedChi2cut_merged"         : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg failed #chi^{2} cut (merged)",    10, [ 0,  10]],

            # "Njets_chi2_lep_failedChi2cut"                : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg failed #chi^{2} cut",             10, [ 0,  10]],
            # "Njets_chi2_lep_failedChi2cut_resolved"       : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg failed #chi^{2} cut (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_lep_failedChi2cut_merged"         : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg failed #chi^{2} cut (merged)",    10, [ 0,  10]],

            # "Njets_chi2_passedMatching"                   : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate",             10, [ 0,  10]],
            # "Njets_chi2_passedMatching_resolved"          : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_passedMatching_merged"            : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate (merged)",    10, [ 0,  10]],

            # "BestChi2_failedMatching"                     : ["#chi^{2} failed matching",                                        50, [ 0, 100]],
            # "BestChi2_failedMatching_resolved"            : ["#chi^{2} failed matching (resolved)",                             50, [ 0, 100]],
            # "BestChi2_failedMatching_merged"              : ["#chi^{2} failed matching (merged)",                               50, [ 0, 100]],

            "Njets_chi2_failedMatching"                   : ["N_{jets} of NotMatched Best #Chi^{2} Candidate",             10, [ 0,  10]],
            "Njets_chi2_failedMatching_resolved"          : ["N_{jets} of NotMatched Best #Chi^{2} Candidate (resolved)",  10, [ 0,  10]],
            "Njets_chi2_failedMatching_merged"            : ["N_{jets} of NotMatched Best #Chi^{2} Candidate (merged)",    10, [ 0,  10]],

            # "BestChi2hadronicleg"                         : ["#chi^{2} of Hadronic leg",                                        50, [ 0, 100]],
            # "BestChi2hadronicleg_resolved"                : ["#chi^{2} of Hadronic leg (resolved)",                             50, [ 0, 100]],
            # "BestChi2hadronicleg_merged"                  : ["#chi^{2} of Hadronic leg (merged)",                               50, [ 0, 100]],

            # "Njets_chi2_had"                              : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg",              10, [ 0,  10]],
            # "Njets_chi2_had_resolved"                     : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg (resolved)",   10, [ 0,  10]],
            # "Njets_chi2_had_merged"                       : ["N_{jets} of Best #Chi^{2} Candidate's hadronic leg (merged)",     10, [ 0,  10]],

            # "BestChi2hadronicleg_passedMatching"          : ["#chi^{2} of Hadronic leg passed matching",                        50, [ 0, 100]],
            # "BestChi2hadronicleg_passedMatching_resolved" : ["#chi^{2} of Hadronic leg passed matching (resolved)",             50, [ 0, 100]],
            # "BestChi2hadronicleg_passedMatching_merged"   : ["#chi^{2} of Hadronic leg passed matching (merged)",               50, [ 0, 100]],

            # "Njets_chi2_had_passedMatching"               : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate's hadronic leg",             10, [ 0,  10]],
            # "Njets_chi2_had_passedMatching_resolved"      : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate's hadronic leg (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_had_passedMatching_merged"        : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate's hadronic leg (merged)",    10, [ 0,  10]],

            # "BestChi2hadronicleg_failedMatching"          : ["#chi^{2} of Hadronic leg failed matching",                        50, [ 0, 100]],
            # "BestChi2hadronicleg_failedMatching_resolved" : ["#chi^{2} of Hadronic leg failed matching (resolved)",             50, [ 0, 100]],
            # "BestChi2hadronicleg_failedMatching_merged"   : ["#chi^{2} of Hadronic leg failed matching (merged)",               50, [ 0, 100]],

            "Njets_chi2_had_failedMatching"               : ["N_{jets} of NotMatched Best #Chi^{2} Candidate's hadronic leg",             10, [ 0,  10]],
            "Njets_chi2_had_failedMatching_resolved"      : ["N_{jets} of NotMatched Best #Chi^{2} Candidate's hadronic leg (resolved)",  10, [ 0,  10]],
            "Njets_chi2_had_failedMatching_merged"        : ["N_{jets} of NotMatched Best #Chi^{2} Candidate's hadronic leg (merged)",    10, [ 0,  10]],

            # "BestChi2leptonicleg"                         : ["#chi^{2} of Leptonic leg",                                        50, [ 0, 100]],
            # "BestChi2leptonicleg_resolved"                : ["#chi^{2} of Leptonic leg (resolved)",                             50, [ 0, 100]],
            # "BestChi2leptonicleg_merged"                  : ["#chi^{2} of Leptonic leg (merged)",                               50, [ 0, 100]],

            # "Njets_chi2_lep"                              : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg",              10, [ 0,  10]],
            # "Njets_chi2_lep_resolved"                     : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg (resolved)",   10, [ 0,  10]],
            # "Njets_chi2_lep_merged"                       : ["N_{jets} of Best #Chi^{2} Candidate's leptonic leg (merged)",     10, [ 0,  10]],

            # "BestChi2leptonicleg_passedMatching"          : ["#chi^{2} of Leptonic leg passed matching",                        50, [ 0, 100]],
            # "BestChi2leptonicleg_passedMatching_resolved" : ["#chi^{2} of Leptonic leg passed matching (resolved)",             50, [ 0, 100]],
            # "BestChi2leptonicleg_passedMatching_merged"   : ["#chi^{2} of Leptonic leg passed matching (merged)",               50, [ 0, 100]],

            # "Njets_chi2_lep_passedMatching"               : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate's leptonic leg",             10, [ 0,  10]],
            # "Njets_chi2_lep_passedMatching_resolved"      : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate's leptonic leg (resolved)",  10, [ 0,  10]],
            # "Njets_chi2_lep_passedMatching_merged"        : ["N_{jets} of CorrectlyMatched Best #Chi^{2} Candidate's leptonic leg (merged)",    10, [ 0,  10]],

            # "BestChi2leptonicleg_failedMatching"          : ["#chi^{2} of Leptonic leg failed matching",                        50, [ 0, 100]],
            # "BestChi2leptonicleg_failedMatching_resolved" : ["#chi^{2} of Leptonic leg failed matching (resolved)",             50, [ 0, 100]],
            # "BestChi2leptonicleg_failedMatching_merged"   : ["#chi^{2} of Leptonic leg failed matching (merged)",               50, [ 0, 100]],

            "Njets_chi2_lep_failedMatching"               : ["N_{jets} of NotMatched Best #Chi^{2} Candidate's leptonic leg",             10, [ 0,  10]],
            "Njets_chi2_lep_failedMatching_resolved"      : ["N_{jets} of NotMatched Best #Chi^{2} Candidate's leptonic leg (resolved)",  10, [ 0,  10]],
            "Njets_chi2_lep_failedMatching_merged"        : ["N_{jets} of NotMatched Best #Chi^{2} Candidate's leptonic leg (merged)",    10, [ 0,  10]],

            # "DrHadB"                                      : ["#DeltaR of hadronic-b",                                           50, [ 0,  10]],
            # "DrHadB_resolved"                             : ["#DeltaR of hadronic-b (resolved)",                                50, [ 0,  10]],
            # "DrHadB_merged"                               : ["#DeltaR of hadronic-b (merged)",                                  50, [ 0,   5]],
            # "DrHadB_passedChi2cut"                        : ["#DeltaR of hadronic-b passed #chi^{2}",                           50, [ 0,  10]],
            # "DrHadB_passedChi2cut_resolved"               : ["#DeltaR of hadronic-b passed #chi^{2} (resolved)",                50, [ 0,  10]],
            # "DrHadB_passedChi2cut_merged"                 : ["#DeltaR of hadronic-b passed #chi^{2} (merged)",                  50, [ 0,   5]],
            # "DrHadB_failedChi2cut"                        : ["#DeltaR of hadronic-b failed #chi^{2}",                           50, [ 0,  10]],
            # "DrHadB_failedChi2cut_resolved"               : ["#DeltaR of hadronic-b failed #chi^{2} (resolved)",                50, [ 0,  10]],
            # "DrHadB_failedChi2cut_merged"                 : ["#DeltaR of hadronic-b failed #chi^{2} (merged)",                  50, [ 0,   5]],
            # "DrLep"                                       : ["#DeltaR of lepton",                                               50, [ 0,  10]],
            # "DrLep_resolved"                              : ["#DeltaR of lepton (resolved)",                                    50, [ 0,  10]],
            # "DrLep_merged"                                : ["#DeltaR of lepton (merged)",                                      50, [ 0, 0.2]],
            # "DrLep_passedChi2cut"                         : ["#DeltaR of lepton passed #chi^{2}",                               50, [ 0,  10]],
            # "DrLep_passedChi2cut_resolved"                : ["#DeltaR of lepton passed #chi^{2} (resolved)",                    50, [ 0,  10]],
            # "DrLep_passedChi2cut_merged"                  : ["#DeltaR of lepton passed #chi^{2} (merged)",                      50, [ 0, 0.2]],
            # "DrLep_failedChi2cut"                         : ["#DeltaR of lepton failed #chi^{2}",                               50, [ 0,  10]],
            # "DrLep_failedChi2cut_resolved"                : ["#DeltaR of lepton failed #chi^{2} (resolved)",                    50, [ 0,  10]],
            # "DrLep_failedChi2cut_merged"                  : ["#DeltaR of lepton failed #chi^{2} (merged)",                      50, [ 0, 0.2]],
    
            # "correctDr"                                   : ["correct_dr",                                                        50, [ 0, 2.0]],
            # "dr_Lep_b"                                    : ["dr contribution of Leptonic b",                                     50, [ 0, 0.4]],
            # "dr_Had_b"                                    : ["dr contribution of Hadronic b",                                     50, [ 0, 0.8]],
            # "dr_Had_q1"                                   : ["dr contribution of Hadronic q1",                                    50, [ 0, 0.8]],
            # "dr_Had_q2"                                   : ["dr contribution of Hadronic q2",                                    50, [ 0, 0.8]],
            # "dr_Lep_nu"                                   : ["dr contribution of Leptonic #nu",                                   50, [ 0, 0.3]],
            # "dr_Lep_l"                                    : ["dr contribution of Leptonic lepton",                                50, [ 0, 0.1]],

            # "correctDr_resolved"                          : ["correct_dr (resolved)",                                             50, [ 0, 2.0]],
            # "dr_Lep_b_resolved"                           : ["dr contribution of Leptonic b (resolved)",                          50, [ 0, 0.4]],
            # "dr_Had_b_resolved"                           : ["dr contribution of Hadronic b (resolved)",                          50, [ 0, 0.4]],
            # "dr_Had_q1_resolved"                          : ["dr contribution of Hadronic q1 (resolved)",                         50, [ 0, 0.4]],
            # "dr_Had_q2_resolved"                          : ["dr contribution of Hadronic q2 (resolved)",                         50, [ 0, 0.4]],
            # "dr_Lep_nu_resolved"                          : ["dr contribution of Leptonic #nu (resolved)",                        50, [ 0, 0.3]],
            # "dr_Lep_l_resolved"                           : ["dr contribution of Leptonic lepton (resolved)",                     50, [ 0, 0.1]],

            # "correctDr_merged"                            : ["correct_dr (merged)",                                               50, [ 0, 2.0]],
            # "dr_Lep_b_merged"                             : ["dr contribution of Leptonic b (merged)",                            50, [ 0, 0.4]],
            # "dr_Had_b_merged"                             : ["dr contribution of Hadronic b (merged)",                            50, [ 0, 0.8]],
            # "dr_Had_q1_merged"                            : ["dr contribution of Hadronic q1 (merged)",                           50, [ 0, 0.8]],
            # "dr_Had_q2_merged"                            : ["dr contribution of Hadronic q2 (merged)",                           50, [ 0, 0.8]],
            # "dr_Lep_nu_merged"                            : ["dr contribution of Leptonic #nu (merged)",                          50, [ 0, 0.3]],
            # "dr_Lep_l_merged"                             : ["dr contribution of Leptonic lepton (merged)",                       50, [ 0, 0.1]],

            # "correctDr_passedChi2cut"                     : ["correct_dr of t#bar{t} passed #chi^{2} cut",                        50, [ 0, 2.0]],
            # "dr_Lep_b_passedChi2cut"                      : ["dr contribution of Leptonic b passed #chi^{2} cut",                 50, [ 0, 0.4]],
            # "dr_Had_b_passedChi2cut"                      : ["dr contribution of Hadronic b passed #chi^{2} cut",                 50, [ 0, 0.8]],
            # "dr_Had_q1_passedChi2cut"                     : ["dr contribution of Hadronic q1 passed #chi^{2} cut",                50, [ 0, 0.8]],
            # "dr_Had_q2_passedChi2cut"                     : ["dr contribution of Hadronic q2 passed #chi^{2} cut",                50, [ 0, 0.8]],
            # "dr_Lep_nu_passedChi2cut"                     : ["dr contribution of Leptonic #nu passed #chi^{2} cut",               50, [ 0, 0.3]],
            # "dr_Lep_l_passedChi2cut"                      : ["dr contribution of Leptonic lepton passed #chi^{2} cut",            50, [ 0, 0.1]],

            # "correctDr_passedChi2cut_resolved"            : ["correct_dr of t#bar{t} passed #chi^{2} cut (resolved)",             50, [ 0, 2.0]],
            # "dr_Lep_b_passedChi2cut_resolved"             : ["dr contribution of Leptonic b passed #chi^{2} cut (resolved)",      50, [ 0, 0.4]],
            # "dr_Had_b_passedChi2cut_resolved"             : ["dr contribution of Hadronic b passed #chi^{2} cut (resolved)",      50, [ 0, 0.4]],
            # "dr_Had_q1_passedChi2cut_resolved"            : ["dr contribution of Hadronic q1 passed #chi^{2} cut (resolved)",     50, [ 0, 0.4]],
            # "dr_Had_q2_passedChi2cut_resolved"            : ["dr contribution of Hadronic q2 passed #chi^{2} cut (resolved)",     50, [ 0, 0.4]],
            # "dr_Lep_nu_passedChi2cut_resolved"            : ["dr contribution of Leptonic #nu passed #chi^{2} cut (resolved)",    50, [ 0, 0.3]],
            # "dr_Lep_l_passedChi2cut_resolved"             : ["dr contribution of Leptonic lepton passed #chi^{2} cut (resolved)", 50, [ 0, 0.1]],

            # "correctDr_passedChi2cut_merged"              : ["correct_dr of t#bar{t} passed #chi^{2} cut (merged)",               50, [ 0, 2.0]],
            # "dr_Lep_b_passedChi2cut_merged"               : ["dr contribution of Leptonic b passed #chi^{2} cut (merged)",        50, [ 0, 0.4]],
            # "dr_Had_b_passedChi2cut_merged"               : ["dr contribution of Hadronic b passed #chi^{2} cut (merged)",        50, [ 0, 0.8]],
            # "dr_Had_q1_passedChi2cut_merged"              : ["dr contribution of Hadronic q1 passed #chi^{2} cut (merged)",       50, [ 0, 0.8]],
            # "dr_Had_q2_passedChi2cut_merged"              : ["dr contribution of Hadronic q2 passed #chi^{2} cut (merged)",       50, [ 0, 0.8]],
            # "dr_Lep_nu_passedChi2cut_merged"              : ["dr contribution of Leptonic #nu passed #chi^{2} cut (merged)",      50, [ 0, 0.3]],
            # "dr_Lep_l_passedChi2cut_merged"               : ["dr contribution of Leptonic lepton passed #chi^{2} cut (merged)",   50, [ 0, 0.1]],

            # "correctDr_failedChi2cut"                     : ["correct_dr of t#bar{t} failed #chi^{2} cut",                        50, [ 0, 2.0]],
            # "dr_Lep_b_failedChi2cut"                      : ["dr contribution of Leptonic b failed #chi^{2} cut",                 50, [ 0, 0.4]],
            # "dr_Had_b_failedChi2cut"                      : ["dr contribution of Hadronic b failed #chi^{2} cut",                 50, [ 0, 0.8]],
            # "dr_Had_q1_failedChi2cut"                     : ["dr contribution of Hadronic q1 failed #chi^{2} cut",                50, [ 0, 0.8]],
            # "dr_Had_q2_failedChi2cut"                     : ["dr contribution of Hadronic q2 failed #chi^{2} cut",                50, [ 0, 0.8]],
            # "dr_Lep_nu_failedChi2cut"                     : ["dr contribution of Leptonic #nu failed #chi^{2} cut",               50, [ 0, 0.5]],
            # "dr_Lep_l_failedChi2cut"                      : ["dr contribution of Leptonic lepton failed #chi^{2} cut",            50, [ 0, 0.1]],

            # "correctDr_failedChi2cut_resolved"            : ["correct_dr of t#bar{t} failed #chi^{2} cut (resolved)",             50, [ 0, 2.0]],
            # "dr_Lep_b_failedChi2cut_resolved"             : ["dr contribution of Leptonic b failed #chi^{2} cut (resolved)",      50, [ 0, 0.4]],
            # "dr_Had_b_failedChi2cut_resolved"             : ["dr contribution of Hadronic b failed #chi^{2} cut (resolved)",      50, [ 0, 0.4]],
            # "dr_Had_q1_failedChi2cut_resolved"            : ["dr contribution of Hadronic q1 failed #chi^{2} cut (resolved)",     50, [ 0, 0.4]],
            # "dr_Had_q2_failedChi2cut_resolved"            : ["dr contribution of Hadronic q2 failed #chi^{2} cut (resolved)",     50, [ 0, 0.4]],
            # "dr_Lep_nu_failedChi2cut_resolved"            : ["dr contribution of Leptonic #nu failed #chi^{2} cut (resolved)",    50, [ 0, 0.3]],
            # "dr_Lep_l_failedChi2cut_resolved"             : ["dr contribution of Leptonic lepton failed #chi^{2} cut (resolved)", 50, [ 0, 0.1]],

            # "correctDr_failedChi2cut_merged"              : ["correct_dr of t#bar{t} failed #chi^{2} cut (merged)",               50, [ 0, 2.0]],
            # "dr_Lep_b_failedChi2cut_merged"               : ["dr contribution of Leptonic b failed #chi^{2} cut (merged)",        50, [ 0, 0.4]],
            # "dr_Had_b_failedChi2cut_merged"               : ["dr contribution of Hadronic b failed #chi^{2} cut (merged)",        50, [ 0, 0.8]],
            # "dr_Had_q1_failedChi2cut_merged"              : ["dr contribution of Hadronic q1 failed #chi^{2} cut (merged)",       50, [ 0, 0.8]],
            # "dr_Had_q2_failedChi2cut_merged"              : ["dr contribution of Hadronic q2 failed #chi^{2} cut (merged)",       50, [ 0, 0.8]],
            # "dr_Lep_nu_failedChi2cut_merged"              : ["dr contribution of Leptonic #nu failed #chi^{2} cut (merged)",      50, [ 0, 0.3]],
            # "dr_Lep_l_failedChi2cut_merged"               : ["dr contribution of Leptonic lepton failed #chi^{2} cut (merged)",   50, [ 0, 0.1]],

            # "DrHadB_matchable"                            : ["#DeltaR of matchable hadronic-b",                                 50, [ 0,  10]],
            # "DrHadB_matchable_resolved"                   : ["#DeltaR of matchable hadronic-b (resolved)",                      50, [ 0,  10]],
            # "DrHadB_matchable_merged"                     : ["#DeltaR of matchable hadronic-b (merged)",                        50, [ 0,   5]],
            # "DrHadB_matchable_passedChi2cut"              : ["#DeltaR of matchable hadronic-b passed #chi^{2} cut",             50, [ 0,  10]],
            # "DrHadB_matchable_passedChi2cut_resolved"     : ["#DeltaR of matchable hadronic-b passed #chi^{2} cut (resolved)",  50, [ 0,  10]],
            # "DrHadB_matchable_passedChi2cut_merged"       : ["#DeltaR of matchable hadronic-b passed #chi^{2} cut (merged)",    50, [ 0,   5]],
            # "DrHadB_matchable_failedChi2cut"              : ["#DeltaR of matchable hadronic-b failed #chi^{2} cut",             50, [ 0,  10]],
            # "DrHadB_matchable_failedChi2cut_resolved"     : ["#DeltaR of matchable hadronic-b failed #chi^{2} cut (resolved)",  50, [ 0,  10]],
            # "DrHadB_matchable_failedChi2cut_merged"       : ["#DeltaR of matchable hadronic-b failed #chi^{2} cut (merged)",    50, [ 0,   5]],
            # "DrLep_matchable"                             : ["#DeltaR of matchable lepton",                                     50, [ 0, 0.2]],
            # "DrLep_matchable_resolved"                    : ["#DeltaR of matchable lepton (resolved)",                          50, [ 0, 0.2]],
            # "DrLep_matchable_merged"                      : ["#DeltaR of matchable lepton (merged)",                            50, [ 0, 0.2]],
            # "DrLep_matchable_passedChi2cut"               : ["#DeltaR of matchable lepton passed #chi^{2} cut",                 50, [ 0, 0.2]],
            # "DrLep_matchable_passedChi2cut_resolved"      : ["#DeltaR of matchable lepton passed #chi^{2} cut (resolved)",      50, [ 0, 0.2]],
            # "DrLep_matchable_passedChi2cut_merged"        : ["#DeltaR of matchable lepton passed #chi^{2} cut (merged)",        50, [ 0, 0.2]],
            # "DrLep_matchable_failedChi2cut"               : ["#DeltaR of matchable lepton failed #chi^{2} cut",                 50, [ 0, 0.2]],
            # "DrLep_matchable_failedChi2cut_resolved"      : ["#DeltaR of matchable lepton failed #chi^{2} cut (resolved)",      50, [ 0, 0.2]],
            # "DrLep_matchable_failedChi2cut_merged"        : ["#DeltaR of matchable lepton failed #chi^{2} cut (merged)",        50, [ 0, 0.2]],
            }

### end User input section i.e. modify each time new set of plots are being generated ------------------------------------------------------------------------------------

### Histogram Information:
# Toggle Stats 
gStyle.SetOptStat("iorme")
# Set canvas width and height in pixels
H = 600
W = 800

# Create canvas
canvas = TCanvas('canvas', 'canvas', W, H)

# Clear text file content
open(txtfile, 'w').close()
print "Cleared contents of", txtfile
print "\n"

# Counter to track progress in output
counter = 1

# Loop through each histogram we want to plot
for hist in hists:
    print "Plotting", hist, "(", counter, "/", len(hists), ")"
    histogram = TH1F("%s"%(hist),"%s"%(hist),hists[hist][1], hists[hist][2][0], hists[hist][2][1])
    histogram.Sumw2()
    # Loop over each sample (different files)
    for sample in sample_names:
        file = TFile("%s/uhh2.AnalysisModuleRunner.MC.%s.root"%(fileDir,sample),"read")
        tree = file.Get("AnalysisTree")
        tree.Draw("%s>>h_%s(%i,%f,%f)"%(hist, hist, hists[hist][1], hists[hist][2][0], hists[hist][2][1]))
        histogram = tree.GetHistogram()

    # Set histogram style
    histogram.SetFillColor(kRed)
    histogram.SetLineColor(0)
    histogram.SetTitle("")
    histogram.GetYaxis().SetTitle("Events")
    histogram.GetXaxis().SetTitle("%s"%(hists[hist][0]).replace("_rebin", ""))
    histogram.Draw("hist")

    # Save histogram file
    if log:
        histogram.SetMinimum(1)
        canvas.SetLogy()
        canvas.SaveAs("%s/%s_log.png"%(plotDirectory, hist))
    else:
        canvas.SaveAs("%s/%s.png"%(plotDirectory, hist))
    canvas.Clear() # Clear canvas for next plot

    ## Begin saving integral and efficiency info to text file --------------------
    if "BestChi2" in hist:
        # Find the bin that corresponds to the cut-value
        cut_value = 0
        cut_value = 30
        cutbin = -1
        for bin in range(0, histogram.GetNbinsX()+1):
            if histogram.GetXaxis().GetBinLowEdge(bin) < cut_value:
                cutbin = bin

        # Compute efficienty info
        Tots = histogram.Integral(1,hists[hist][1]+1)
        Npassed = histogram.Integral(1, cutbin)
        Efficiency = Npassed/Tots

        print "Npassed =", Npassed, "and Tots =", Tots
        print "Efficiency =", Efficiency

        # Save Integral info to text file
        with open(txtfile, 'a') as f:
            f.write('%s:\n'%(hists[hist][0]))
            f.write('Npassed (bins 1 to %i): %i\n'%(cutbin, int(Npassed)) )
            f.write('Total (bins 1 to %i): %i\n'%(hists[hist][1]+1, int(Tots)) )
            f.write('Efficiency is %f \n'%(Efficiency) )
            print "Appended Integral and efficiency info to", txtfile
            f.write('\n')
            f.close()
    else:
        Tots = histogram.Integral(1,hists[hist][1]+1)
        print "Total (bins 1 to %i): %i\n"%(hists[hist][1]+1, int(Tots))
        with open(txtfile, 'a') as f:
            f.write('%s:\n'%(hists[hist][0]))
            f.write('Total (bins 1 to %i): %i\n'%(hists[hist][1]+1, int(Tots)) )
            f.write('\n')
            f.close()
    ## End saving integral and efficiency info to text file --------------------
            
    counter += 1
    print "\n" # Seperates each histogram with a new line in the terminal output

print "Hope that worked, good luck with the slides!"