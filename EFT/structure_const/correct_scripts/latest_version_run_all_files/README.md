# Description of relevant scripts

1. calc_structure_constants.py

   - Now processes all ROOT files in a directory
   - Uses parallel processing for faster execution
   - Saves structure constants in a "structure_constants" folder
   - Each output file is named after the input file

2. plot_deltaY_shape.py *and* plot_deltaPhi_sigmaPhi.py:

   - Added option to save histograms to ROOT files
   - Can now be used to generate both plots and histogram files

3. process_all_files.py

   - Finds all ROOT files and matches them with their structure constants files
   - Runs the plotting scripts for each pair
   - Saves individual histograms in separate directories
   - Combines all histograms using hadd

4. create_combined_plots.py

   - This script creates combined ratio plots from existing ROOT files.
   - Reads the individual ratio histograms from ROOT files (output of process_all_files.py) and creates
   - combined plots with all 16 Wilson coefficients on the same canvas.

5. extract_plots_to_pdf.py

   - This script extracts the ratio plots from the ROOT files (output of process_all_files.py) and saves them as PDFs.

--------------------------------

## Usage

**python calc_structure_constants.py** /path/to/your/root/files/ 8

/data/dust/group/cms/zprime-uhh/Analysis_EFT_UL17/muon/workdir_Analysis_EFT_UL17_muon_semilepton_deltayreco_deltaPhi_sigmaPhi_ttree

python **process_all_files.py** /path/to/your/root/files/ structure_constants_all_files

python **process_all_files.py** /data/dust/group/cms/zprime-uhh/Analysis_EFT_UL17/muon/workdir_Analysis_EFT_UL17_muon_semilepton_deltayreco_deltaPhi_sigmaPhi_ttree/ ../../correct_scripts/structure_constants_all_files_2017/

This will:

- Match each ROOT file with its corresponding structure constants file
- Run the plotting scripts for each pair
- Save individual histograms in the "combined_output" directory
- Combine all histograms into three files:

combined_deltaY.root
combined_deltaPhi.root
combined_sigmaPhi.root

Optional flags:
--custom-wc: Use custom Wilson coefficient values defined in the scripts
--delta-phi-only: Only process Delta Phi
--sigma-phi-only: Only process Sigma Phi
--delta-y-only: Only process Delta Y
--output-dir DIR: Specify a custom output directory (default: "combined_output")

For each input ROOT file (example_file.root) with its matching structure constants file, the script will create:
Delta Y histograms:
Saved as: combined_output/deltaY/example_file_deltaY.root
Delta Phi histograms:
Saved as: combined_output/deltaPhi/example_file_deltaPhi.root
Sigma Phi histograms:
Saved as: combined_output/sigmaPhi/example_file_sigmaPhi.root

combined_output/plots/example_file/deltaY_plots/PDFs
combined_output/plots/example_file/deltaPhi_plots/PDFs
combined_output/plots/example_file/sigmaPhi_plots/PDFs

--------------------------------

For each input file, the scripts save:

- Individual Plots in ROOT File

For each Wilson coefficient (WC), the individual plots are saved as:

1) SM Histograms:
hDeltaY_SM_<wc_name> - Standard Model distribution for Delta Y
hDeltaPhi_SM_<wc_name> - Standard Model distribution for Delta Phi
hSigmaPhi_SM_<wc_name> - Standard Model distribution for Sigma Phi

1) EFT Histograms:
hDeltaY_EFT_<wc_name> - EFT distribution for Delta Y
hDeltaPhi_EFT_<wc_name> - EFT distribution for Delta Phi
hSigmaPhi_EFT_<wc_name> - EFT distribution for Sigma Phi

1) Normalized Histograms:
hDeltaY_SM_norm_<wc_name> - Normalized SM distribution
hDeltaY_EFT_norm_<wc_name> - Normalized EFT distribution
(Similar for Delta Phi and Sigma Phi)

1) Individual Ratio Histograms:
hRatio_norm_deltaY_<wc_name> - Ratio of EFT/SM for Delta Y
hRatio_norm_deltaPhi_<wc_name> - Ratio of EFT/SM for Delta Phi
hRatio_norm_sigmaPhi_<wc_name> - Ratio of EFT/SM for Sigma Phi
These individual histograms preserve error bars in the ROOT file.

- To create combined ratio plots, run:
python create_combined_plots.py combined_output/root_files/ combined_output/plots/
single root file: python create_combined_plots.py combined_output/root_files/uhh2.AnalysisModuleRunner.MC.MC_EFT_Mttbar_0-700_UL17_8_plots.root --output-dir combined_plots/Mttbar_0-700
    directory: python create_combined_plots.py combined_output/root_files/ --output-dir combined_plots

- To create individual plots, run:
Single root file:
python plot_deltaY_shape.py combined_output/root_files/ combined_output/plots/
python plot_deltaPhi_sigmaPhi.py combined_output/root_files/ combined_output/plots/

OR

hadd root files:
python extract_plots_to_pdf.py combined_output/root_files/ combined_output/plots/
python extract_plots_to_pdf.py combined_output_v1/root_files/uhh2.AnalysisModuleRunner.MC.MC_EFT_Mttbar_0-700_UL17_36_all_plots.root --output-dir plots_with_extractplotscript
