from ROOT import *
import os
import sys
import math
from optparse import OptionParser
from ROOT import TH1D
import math
      

parser = OptionParser()
parser.add_option("-y", "--year", dest="year", help="Specify the year (UL18, UL17, preUL16, postUL16)", type='str')
parser.add_option("-m", "--mass_range", dest="mass_range", help="Specify the mass range (0_500, 500_750, 750-1000, 1000-1500, 1500Inf)", type='str')
parser.add_option("-l", "--lepton_flavor", dest="lepton_flavor", help="Specify the lepton flavor (ele, muon)", type='str')
parser.add_option("-r", "--region", dest="region", help="Specify the region (SR, CR1, CR2)", type='str')
parser.add_option("-s", "--region_score", dest="region_score", help="Specify the region score (0, 1, 2)", type='str')


(options, args) = parser.parse_args()

year = options.year if options.year else "UL17" 
mass_range = options.mass_range if options.mass_range else "0_500" 
lepton_flavor = options.lepton_flavor if options.lepton_flavor else "muon" 
# region = options.region if options.region else "SR"
# region = options.region if options.region else "CR1"
region = options.region if options.region else "CR2"

region_score = options.region_score if options.region_score else "0"



# finalState = options.channel

inputDir = "/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_DNN/{}/{}/nominal/".format(year, lepton_flavor)
combine_file_name = 'dY_{}_{}_{}_{}.root'.format(year, lepton_flavor, mass_range, region)
combine_file = TFile(combine_file_name, 'RECREATE')
# stackList = {"TTbar", "WJets", "DY", "ST", "Diboson", "QCD", "DATA"}
stackList = ["TTbar", "ST", "Others", "DATA"]


############ 
# ----------------- other systematics except PDF, JER/JEC, murmuf ----------------
############ 

debug = False

systematic_name_mapping = {
    "mu_reco": "muonReco",
    "pu": "pu",
    "prefiring": "prefiringWeight",
    "mu_id_stat" : "muonIDStat",
    "mu_id_syst" : "muonIDSyst",
    "mu_iso_stat" : "muonIsoStat",
    "mu_iso_syst" : "muonIsoSyst",
    "mu_trigger_stat" : "muonTriggerStat",
    "mu_trigger_syst" : "muonTriggerSyst",
    "ele_id" : "eleID", 
    "ele_trigger": "eleTrigger",
    "ele_reco": "eleReco",
    "isr": "isr", 
    "fsr": "fsr", 
    "btag_cferr1": "btagCferr1", 
    "btag_cferr2": "btagCferr2", 
    "btag_hf": "btagHf",  
    "btag_hfstats1": "btagHfstats1",
    "btag_hfstats2": "btagHfstats2", 
    "btag_lf": "btagLf", 
    "btag_lfstats1": "btagLfstats1",
    "btag_lfstats2": "btagLfstats2", 
    "ttag_corr": "ttagCorr", 
    "ttag_uncorr": "ttagUncorr", 
    "tmistag": "tmistag"  
    
}

# Define lepton-specific systematics
muon_systematics = ["mu_reco", "mu_id_stat", "mu_id_syst", "mu_iso_stat", "mu_iso_syst", "mu_trigger_stat", "mu_trigger_syst"]
electron_systematics = ["ele_id", "ele_trigger", "ele_reco"]

# reads nominal histograms and computes projections for the "TTbar" sample. 
# for other samples, it reads nominal histograms and writes them to the output file. 
# It also processes and writes systematic variations, using systematic_name_mapping

if (debug): print(" ----------------- all systematics except PDF, JER/JEC, murmuf ----------------")

for sample in stackList:
    inFile = TFile.Open(inputDir + "{}.root".format(sample), "READ")
    if not inFile:
        print("Input file for {} not found.".format(sample))
        continue
    
    combine_file.cd()

    
    if sample == "DATA":
        # print("in data loop, sample: ", sample)
        if (debug): print("Processing sample for other sys: ", sample)
        data_obs = inFile.Get("DeltaY_reco_{}_{}_General/DeltaY_reco".format(mass_range,region)).Clone("data_obs")
        data_obs.Write("data_obs")
        
    elif sample == "TTbar":
        if (debug): print("Processing sample for other sys: ", sample)
        path_nom = "DeltaY_reco_SystVariations_{}_{}/DeltaY_tt".format(mass_range, region)
        # Retrieving the response matrix directly
        Matrix = inFile.Get(path_nom)
        if not Matrix:
            print("Response matrix not found for TTbar.")
            continue
        
        ProjX_1 = Matrix.ProjectionX("TTbar_1", 1, 1)
        ProjX_2 = Matrix.ProjectionX("TTbar_2", 2, 2)

        ProjX_1.GetXaxis().SetTitle("#Delta_Y_{reco}")
        ProjX_2.GetXaxis().SetTitle("#Delta_Y_{reco}")

        nominal_projections = [ProjX_1, ProjX_2]
        
        Matrix.Write()
        ProjX_1.Write("TTbar_1")
        ProjX_2.Write("TTbar_2")
        
        
        for sys, new_sys_name in systematic_name_mapping.items():
            # Skip electron systematics for muon flavor
            if lepton_flavor == "muon" and sys in electron_systematics:
                if (debug): print("Skipping electron systematic {} for muon flavor".format(sys))
                continue
            # Skip muon systematics for electron flavor
            if lepton_flavor == "ele" and sys in muon_systematics:
                if (debug): print("Skipping muon systematic {} for electron flavor".format(sys))
                continue
                
            if (debug): print(sys)
            for variation in ["up", "down"]:
                
                # Retrieve matrices for systematic variations
                matrix_path = "DeltaY_reco_SystVariations_{}_{}/DeltaY_{}_{}_tt".format(mass_range, region, sys,variation)
                Matrix_var = inFile.Get(matrix_path)
                if not Matrix_var:
                    print("Response matrix for systematic {} variation {} not found.").format(sys, variation)
                    continue
                
                output_hist_1 = "TTbar_1_{}{}".format(new_sys_name, variation.capitalize())
                output_hist_2 = "TTbar_2_{}{}".format(new_sys_name, variation.capitalize())
                
                ProjX_1_sys = Matrix_var.ProjectionX(output_hist_1, 1, 1)
                ProjX_2_sys = Matrix_var.ProjectionX(output_hist_2, 2, 2)

                ProjX_1_sys.GetXaxis().SetTitle("#Delta_Y_{reco}")
                ProjX_2_sys.GetXaxis().SetTitle("#Delta_Y_{reco}")

                ProjX_1_sys.Write(output_hist_1)
                ProjX_2_sys.Write(output_hist_2)
                
                # print("Nominal Projection mean:", ProjX_1.GetBinContent(1))
                # print("Systematic Projection mean:", ProjX_1_sys.GetBinContent(1))

    else:
        if (debug): print("Processing sample for other sys: ", sample)
        
        h_nominal = inFile.Get("DeltaY_reco_SystVariations_{}_{}/DeltaY".format(mass_range, region)).Clone(sample)
        h_nominal.Write(sample)
        
        
        for sys, new_sys_name in systematic_name_mapping.items():
            # Skip electron systematics for muon flavor
            if lepton_flavor == "muon" and sys in electron_systematics:
                if (debug): print("Skipping electron systematic {} for muon flavor".format(sys))
                continue
            # Skip muon systematics for electron flavor
            if lepton_flavor == "ele" and sys in muon_systematics:
                if (debug): print("Skipping muon systematic {} for electron flavor".format(sys))
                continue
                
            if (debug): print(sys)
            
            for variation in ["up", "down"]:
                sys_hist_name = "DeltaY_{}_{}".format(sys, variation)
                sys_hist = inFile.Get("DeltaY_reco_SystVariations_{}_{}/".format(mass_range, region) + sys_hist_name)
                if sys_hist:
                    output_hist_name = "{}_{}{}".format(sample, new_sys_name, variation.capitalize())
                    sys_hist.Clone(output_hist_name).Write()
            
         
    inFile.Close()

# processes each sample in v_samples to create envelope histograms for each variation in v_variations. 
# histograms are normalized based on the first bin content and written to the DeltaY directory in the output ROOT file

# calculates normalization factors, and then computes the maximum and minimum variations for each bin. The resulting histograms are named sample_murmufUp and sample_murmufDown

############ 
# ----------------- murmuf ----------------

# 1) Retrieve histograms from the input file: nominal, upup, upnone
# 2) Calculate normalization scales for each variation by dividing the first bin content of each variation histogram by the first bin content of the nominal histogram.
# 3) Apply these normalization scales to the corresponding histograms from a vector of murmuf variations (like v_murmuf_upup, v_murmuf_upnone, etc.).
# 4) For each set of variation histograms, find the maximum and minimum values in each bin among all variations, create two new histograms: hist_scaleUp and hist_scaleDown.
# 5) Name these histograms depending on the specific convention
############ 

def calculate_envelope(projection, nominal_projection, norm, up_name, down_name):
    hist_scale_up = nominal_projection.Clone(up_name)
    hist_scale_down = nominal_projection.Clone(down_name)
    
    for bin_idx in range(1, nominal_projection.GetNbinsX() + 1):
        scaled_val = projection.GetBinContent(bin_idx) / norm
        max_val = max(nominal_projection.GetBinContent(bin_idx), scaled_val)
        min_val = min(nominal_projection.GetBinContent(bin_idx), scaled_val)
        hist_scale_up.SetBinContent(bin_idx, max_val)
        hist_scale_down.SetBinContent(bin_idx, min_val)

    return hist_scale_up, hist_scale_down

def getEnvelope(inputDir, v_samples_ttbar, v_variations, combine_file):
    
    if (debug): print(" ----------------- murmuf processing ----------------")

    
    for sample in v_samples_ttbar:
        inFile = TFile.Open(inputDir + "{}.root".format(sample), "READ")       
        if not inFile:
            print("Input file for {} not found.".format(sample))
            continue
        combine_file.cd()
        
        
        if sample == "TTbar":
        
            if (debug): print("Processing sample for murmuf: ", sample)
            
            scaled_histograms_up_1 = {}
            scaled_histograms_down_1 = {}
            scaled_histograms_up_2 = {}
            scaled_histograms_down_2 = {}
            
            hist_scaleUp_1 = ProjX_1.Clone("TTbar_1_murmufUp")
            hist_scaleDown_1 = ProjX_1.Clone("TTbar_1_murmufDown")
            hist_scaleUp_2 = ProjX_2.Clone("TTbar_2_murmufUp")
            hist_scaleDown_2 = ProjX_2.Clone("TTbar_2_murmufDown")
        
            for variation in v_variations:
                path_var = "DeltaY_reco_SystVariations_{}_{}/DeltaY_murmuf_{}_tt".format(mass_range, region,variation)
                matrix_var = inFile.Get(path_var)
               
                if not matrix_var:
                    print("Missing matrix for variation: {}".format(variation))
                
                projection_1 = matrix_var.ProjectionX("TTbar_1_murmuf_{}".format(variation), 1, 1)
                projection_2 = matrix_var.ProjectionX("TTbar_2_murmuf_{}".format(variation), 2, 2)
                
                norm_1 = projection_1.GetBinContent(1) / ProjX_1.GetBinContent(1)
                norm_2 = projection_2.GetBinContent(1) / ProjX_2.GetBinContent(1)
                
                scaled_histograms_up_1[variation] = projection_1.Clone("TTbar_1_{}_scaled_up".format(variation))
                scaled_histograms_down_1[variation] = projection_1.Clone("TTbar_1_{}_scaled_down".format(variation))
                scaled_histograms_up_2[variation] = projection_2.Clone("TTbar_2_{}_scaled_up".format(variation))
                scaled_histograms_down_2[variation] = projection_2.Clone("TTbar_2_{}_scaled_down".format(variation))
                
                scaled_histograms_up_1[variation].Scale(1 / norm_1)
                scaled_histograms_down_1[variation].Scale(1 / norm_1)
                scaled_histograms_up_2[variation].Scale(1 / norm_2)
                scaled_histograms_down_2[variation].Scale(1 / norm_2)
                
                # print("Nominal Projection mean:", ProjX_1.GetMean())
                # print("Systematic Projection mean:", projection_1.GetMean())


            for bin_idx in range(1, ProjX_1.GetNbinsX() + 1):
                max_val_1 = ProjX_1.GetBinContent(bin_idx)
                min_val_1 = ProjX_1.GetBinContent(bin_idx)
                max_val_2 = ProjX_2.GetBinContent(bin_idx)
                min_val_2 = ProjX_2.GetBinContent(bin_idx)

                for key in scaled_histograms_up_1:
                    max_val_1 = max(max_val_1, scaled_histograms_up_1[key].GetBinContent(bin_idx))
                for key in scaled_histograms_down_1:    
                    min_val_1 = min(min_val_1, scaled_histograms_down_1[key].GetBinContent(bin_idx))

                for key in scaled_histograms_up_2:
                    max_val_2 = max(max_val_2, scaled_histograms_up_2[key].GetBinContent(bin_idx))
                for key in scaled_histograms_down_2:    
                    min_val_2 = min(min_val_2, scaled_histograms_down_2[key].GetBinContent(bin_idx))

                hist_scaleUp_1.SetBinContent(bin_idx, max_val_1)
                hist_scaleDown_1.SetBinContent(bin_idx, min_val_1)
                hist_scaleUp_2.SetBinContent(bin_idx, max_val_2)
                hist_scaleDown_2.SetBinContent(bin_idx, min_val_2)

            hist_scaleUp_1.Write()
            hist_scaleDown_1.Write()
            hist_scaleUp_2.Write()
            hist_scaleDown_2.Write()
            
            


        inFile.Close()


############ 
# ----------------- PDFs (100) ----------------
############ 

# Function to calculate and write RMS based histograms
def calculate_and_write_rms_histograms(projections, nominal_projection, hist_name):
    
    pdf_up_name = "{}_pdfUp".format(hist_name)
    pdf_down_name = "{}_pdfDown".format(hist_name)
    
    hist_pdf_up = TH1D(pdf_up_name, "", nominal_projection.GetNbinsX(), nominal_projection.GetXaxis().GetXmin(), nominal_projection.GetXaxis().GetXmax())
    hist_pdf_down = TH1D(pdf_down_name, "", nominal_projection.GetNbinsX(), nominal_projection.GetXaxis().GetXmin(), nominal_projection.GetXaxis().GetXmax())

    for bin_idx in range(1, nominal_projection.GetNbinsX() + 1):
        nominal_val = nominal_projection.GetBinContent(bin_idx)
        values = projections.GetBinContent(bin_idx)
        rms = math.sqrt(sum((x - nominal_val) ** 2 for x in values) / len(values))
        
        hist_pdf_up.SetBinContent(bin_idx, nominal_val + rms)
        hist_pdf_down.SetBinContent(bin_idx, nominal_val - rms)

    combine_file.cd()
    hist_pdf_up.Write()
    hist_pdf_down.Write()


def processPDF(inputDir, v_samples_ttbar, combine_file):
    if (debug): print(" ----------------- PDFs processing ----------------")
    
    
    for sample in v_samples_ttbar:
        inFile = TFile.Open(inputDir + "{}.root".format(sample), "READ")
        if not inFile:
            print("Input file for {} not found.".format(sample))
            continue
        combine_file.cd()


        if sample == "TTbar":
            if (debug): print("Processing TTbar for PDF: ", sample)
            
            # Step 1: Retrieve the nominal matrix and create projections
            
            # Process PDF variations
            pdf_projections_1 = []
            pdf_projections_2 = []
            
            for i in range(1, 101):  
                weight_matrix_path = "DeltaY_reco_PDFVariations_{}_{}/DeltaY_PDF_RM_{}".format(mass_range,region, i)
                weight_matrix = inFile.Get(weight_matrix_path)
                if not weight_matrix:
                    print("PDF weight matrix {} not found.").format(i)
                    continue
                
                
                pdf_projections_1.append(weight_matrix.ProjectionX("TTbar_1_pdf_{}".format(i), 1, 1))
                pdf_projections_2.append(weight_matrix.ProjectionX("TTbar_2_pdf_{}".format(i), 2, 2))
                
            # Calculate RMS and write histograms for both projections
            # Step 2: normalization & rms calculation
            norm_scales_1 = []
            norm_scales_2 = []
                
            for i in range(100):
                norm_scale_1 = pdf_projections_1[i].GetBinContent(1) / ProjX_1.GetBinContent(1)
                norm_scales_1.append(norm_scale_1)
                pdf_projections_1[i].Scale(1. / norm_scale_1)
                
                norm_scale_2 = pdf_projections_2[i].GetBinContent(1) / ProjX_2.GetBinContent(1)
                norm_scales_2.append(norm_scale_2)
                pdf_projections_2[i].Scale(1. / norm_scale_2)
                
            hist_pdfUp_1 = TH1D("TTbar_1_pdfUp", "", ProjX_1.GetNbinsX(), ProjX_1.GetXaxis().GetXmin(), ProjX_1.GetXaxis().GetXmax())
            hist_pdfDown_1 = TH1D("TTbar_1_pdfDown", "", ProjX_1.GetNbinsX(), ProjX_1.GetXaxis().GetXmin(), ProjX_1.GetXaxis().GetXmax())
            hist_pdfUp_2 = TH1D("TTbar_2_pdfUp", "", ProjX_2.GetNbinsX(), ProjX_2.GetXaxis().GetXmin(), ProjX_2.GetXaxis().GetXmax())
            hist_pdfDown_2 = TH1D("TTbar_2_pdfDown", "", ProjX_2.GetNbinsX(), ProjX_2.GetXaxis().GetXmin(), ProjX_2.GetXaxis().GetXmax())


            for bin_idx in range(1, ProjX_1.GetNbinsX() + 1):
                sum_bins_1 = 0.
                for i in range(100):
                    bin_content_1_pdf = pdf_projections_1[i].GetBinContent(bin_idx)
                    sum_bins_1 += (bin_content_1_pdf - ProjX_1.GetBinContent(bin_idx)) ** 2
                rms_1 = math.sqrt(sum_bins_1 / 100)
                
        
                value_up_1 = ProjX_1.GetBinContent(bin_idx) + rms_1
                value_down_1 = ProjX_1.GetBinContent(bin_idx) - rms_1
                
                hist_pdfUp_1.SetBinContent(bin_idx, value_up_1 )
                hist_pdfDown_1.SetBinContent(bin_idx, value_down_1)
            
            for bin_idx in range(1, ProjX_1.GetNbinsX() + 1):
                sum_bins_2 = 0.
                for i in range(100):
                    bin_content_2_pdf = pdf_projections_2[i].GetBinContent(bin_idx)
                    sum_bins_2 += (bin_content_2_pdf - ProjX_2.GetBinContent(bin_idx)) ** 2
                
                rms_2 = math.sqrt(sum_bins_2 / 100)
                
                value_up_2 = ProjX_2.GetBinContent(bin_idx) + rms_2
                value_down_2 = ProjX_2.GetBinContent(bin_idx) - rms_2
                
                hist_pdfUp_2.SetBinContent(bin_idx, value_up_2)
                hist_pdfDown_2.SetBinContent(bin_idx, value_down_2)
            
            hist_pdfUp_1.Write()
            hist_pdfDown_1.Write()
            hist_pdfUp_2.Write()
            hist_pdfDown_2.Write()

                
        inFile.Close()
            
            
# Diboson Sample:
# For Diboson, I am cloning the nominal histogram for each of the 100 PDF variations and naming them Up or Down variation. So, this sample is without specific PDF variations.

# TTbar Sample:
# I am creating a matrix from four quadrants, calculating the projection for each, and then applying each of the 100 PDF variations.
# For each PDF variation, I create a matrix, calculate its projections, and then calculate the deviation from the nominal projection for both projections (TTbar_1 and TTbar_2).

# Other Samples:
# For samples other than Diboson and TTbar, I am looping through each of the 100 PDF variations and calculate the Up and Down variations based on the deviation from the nominal histogram.

############ 
# ----------------- JER/JEC ----------------
############ 


def processJERJEC(inputDir, v_samples_ttbar, combine_file, sys_variations, year):
    if (debug): print(" ----------------- JER/JEC processing for year {} ----------------".format(year))

    for sys_variation in sys_variations: 
        if (debug): print(" --- processing {} variation for year {} ---".format(sys_variation, year))
        for sample in v_samples_ttbar:
            sys_file = TFile.Open("/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_DNN/{}/{}/workdir_AnalysisDNN_{}_{}_{}/{}.root".format(year, lepton_flavor, year, lepton_flavor, sys_variation, sample), "READ")
            if not sys_file:
                if (debug): print("Input file for {} variation {} not found.".format(sample, sys_variation))
                continue
            
            combine_file.cd()

            if sample == "TTbar":
                if (debug): print("Processing TTbar for {} in year {}".format(sys_variation, year))
                
                matrix_path_jerjec = "DeltaY_reco_{}_{}_General/response_matrix".format(mass_range, region)
                Matrix = sys_file.Get(matrix_path_jerjec)
                if not Matrix:
                    print("Matrix not found for TTbar for year {}.".format(year))
                else:
                    if "JER_up" in sys_variation or "JER_down" in sys_variation:
                        output_name_1 = "TTbar_1_{}_{}".format(sys_variation.split('_')[0].lower(), year + sys_variation.split('_')[1].capitalize())
                        output_name_2 = "TTbar_2_{}_{}".format(sys_variation.split('_')[0].lower(), year + sys_variation.split('_')[1].capitalize())

                    elif "JEC_up" in sys_variation or "JEC_down" in sys_variation:
                        output_name_1 = "TTbar_1_{}".format(sys_variation.split('_')[0].lower() + sys_variation.split('_')[1].capitalize())
                        output_name_2 = "TTbar_2_{}".format(sys_variation.split('_')[0].lower() + sys_variation.split('_')[1].capitalize())
                    
                    
                    # Create unique output names for each year
                    
                    # Generate the JER/JEC projections for each year
                    ProjX_1_jerjec = Matrix.ProjectionX(output_name_1, 1, 1)
                    ProjX_2_jerjec = Matrix.ProjectionX(output_name_2, 2, 2)
                    
                    # Write the year-specific JER projections to the combine file
                    ProjX_1_jerjec.Write(output_name_1)
                    ProjX_2_jerjec.Write(output_name_2)

#to process individual JEC sources
def processJECSources(inputDir, v_samples_ttbar, combine_file, jec_sources, year):
    if (debug): print(" ----------------- JEC Sources processing for year {} ----------------".format(year))

    for jec_source in jec_sources:
        for direction in ["up", "down"]:
            sys_variation = "JEC_{}_{}".format(jec_source, direction)
            if (debug): print(" --- processing {} variation for year {} ---".format(sys_variation, year))
            
            for sample in v_samples_ttbar:
                sys_file = TFile.Open("/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_JEC/{}/{}/workdir_AnalysisDNN_{}_{}_JEC_{}_{}/{}.root".format(
                                                                                                                                    year, lepton_flavor, year, lepton_flavor, direction, jec_source, sample), "READ")
                
                if not sys_file:
                    if (debug): print("Input file for JEC source {} direction {} not found.".format(jec_source, direction))
                    continue
                
                combine_file.cd()

                if sample == "TTbar":
                    if (debug): print("Processing TTbar for JEC source {} direction {} in year {}".format(jec_source, direction, year))
                    
                    matrix_path_jec = "DeltaY_reco_{}_{}_General/response_matrix".format(mass_range, region)
                    Matrix = sys_file.Get(matrix_path_jec)
                    if not Matrix:
                        print("Matrix not found for TTbar JEC source {} direction {} for year {}.".format(jec_source, direction, year))
                    else:
                        # Format JEC source name for histogram naming
                        formatted_source = jec_source
                        output_name_1 = "TTbar_1_jec{}{}".format(formatted_source, direction.capitalize())
                        output_name_2 = "TTbar_2_jec{}{}".format(formatted_source, direction.capitalize())
                        
                        # Generate the JEC source projections
                        ProjX_1_jec = Matrix.ProjectionX(output_name_1, 1, 1)
                        ProjX_2_jec = Matrix.ProjectionX(output_name_2, 2, 2)
                        
                        # Write the JEC source projections to the combine file
                        ProjX_1_jec.Write(output_name_1)
                        ProjX_2_jec.Write(output_name_2)

############ 
# ----------------- hdamp ----------------
############ 

def processhdamp(inputDir, v_samples_ttbar, combine_file, sys_variations_hdamp, year):
    if (debug): print(" ----------------- hdamp processing for year {} ----------------".format(year))

    for sys_variation in sys_variations_hdamp: 
        if (debug): print(" --- processing {} variation for year {} ---".format(sys_variation, year))
        for sample in v_samples_ttbar:
            sys_file = TFile.Open("/data/dust/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/output_DNN/{}/{}/workdir_AnalysisDNN_{}_{}_{}/{}.root".format(year, lepton_flavor, year, lepton_flavor, sys_variation, sample), "READ")
            # print("file path: ", sys_file)
            if not sys_file:
                if (debug): print("Input file for {} variation {} not found.".format(sample, sys_variation))
                continue
            
            combine_file.cd()

            if sample == "TTbar":
                if (debug): print("Processing TTbar for {} in year {}".format(sys_variation, year))
                
                matrix_path_hdamp = "DeltaY_reco_{}_{}_General/response_matrix".format(mass_range, region)
                Matrix = sys_file.Get(matrix_path_hdamp)
                if not Matrix:
                    print("Matrix not found for TTbar for year {}.".format(year))
                else:
                    output_name_1 = "TTbar_1_{}".format(sys_variation.split('_')[0].lower() + sys_variation.split('_')[1].capitalize())
                    output_name_2 = "TTbar_2_{}".format(sys_variation.split('_')[0].lower() + sys_variation.split('_')[1].capitalize())
                    
                    
                    # Create unique output names for each year
                    
                    # Generate the hdamp projections for each year
                    ProjX_1_hdamp = Matrix.ProjectionX(output_name_1, 1, 1)
                    ProjX_2_hdamp = Matrix.ProjectionX(output_name_2, 2, 2)
                    
                    # Write the year-specific hdamp projections to the combine file
                    ProjX_1_hdamp.Write(output_name_1)
                    ProjX_2_hdamp.Write(output_name_2)

# JEC sources as provided
jec_sources = [
    "AbsoluteStat",
    "AbsoluteScale",
    "AbsoluteMPFBias",
    "FlavorQCD",
    "Fragmentation",
    "PileUpDataMC",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF",
    "PileUpPtRef",
    "RelativeFSR",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeBal",
    "RelativeSample",
    "RelativeStatEC",
    "RelativeStatFSR",
    "RelativeStatHF",
    "SinglePionECAL",
    "SinglePionHCAL",
    "TimePtEta"
]

sys_variations = ["JEC_up", "JEC_down", "JER_up", "JER_down"]
sys_variations_hdamp = ["hdamp_up", "hdamp_down"]

# v_samples = ["TTbar", "WJets", "ST", "QCD", "DY", "Diboson"]
v_samples = ["TTbar", "ST", "Others"]
v_samples_ttbar = ["TTbar"]

v_variations = ["upup", "upnone", "noneup", "nonedown", "downnone", "downdown"]

getEnvelope(inputDir, v_samples_ttbar, v_variations, combine_file)

processPDF(inputDir, v_samples_ttbar, combine_file)

# Call the JER/JEC processing for each year
processJERJEC(inputDir, v_samples_ttbar, combine_file, sys_variations, year)

# Call the new JEC sources processing
processJECSources(inputDir, v_samples_ttbar, combine_file, jec_sources, year)

# Call the hdamp processing 
# processhdamp(inputDir, v_samples_ttbar, combine_file, sys_variations_hdamp, year)


combine_file.Close()



#  For TTbar samples, the nominal histograms are not directly in the input file. The projection histograms I making in the script will be my 2 nominal histograms (TTbar_1 and TTbar_2). I want to save these histograms in the output file as nominal histograms. Additionally, I would like to apply the w