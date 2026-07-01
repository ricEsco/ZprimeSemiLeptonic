# -*- coding: utf-8 -*-
import ROOT

# Store histogram names and corresponding integer multipliers
coeff_multipliers = {
    "cos_theta1k_antiLep": 2,
    "cos_theta1r_antiLep": 2,
    "cos_theta1n_antiLep": 2,
    "cos_theta2k_Lep": -2,
    "cos_theta2r_Lep": -2,
    "cos_theta2n_Lep": -2,
    "cos_theta1k": 5,
    "cos_theta1r": 5,
    "cos_theta1n": 5,
    "cos_theta2k": 5,
    "cos_theta2r": 5,
    "cos_theta2n": 5,
    "Cnn": -10,
    "Cnr": -10,
    "Cnk": -10,
    "Crn": -10,
    "Crr": -10,
    "Crk": -10,
    "Ckn": -10,
    "Ckr": -10,
    "Ckk": -10,
    "cHel":                              5,
    "cHel_Mtt300_400":                   5,
    "cHel_Mtt300_400_betaLT0p9":         5,
    "cHel_P3n":                          5,
    "cHel_P3n_Mtt800_Inf":               5,
    "cHel_P3n_Mtt800_Inf_cosThetaLT0p4": 5,
}

def compute_fb_asymmetry(histogram, multiplier):
    
    # Ensure the histogram has an even number of bins
    num_bins = histogram.GetNbinsX()
    if num_bins % 2 != 0:
        raise ValueError("Histogram must have an even number of bins.")
    half = num_bins // 2

    # Initialize sums for forward and backward hemispheres
    sum_F = 0
    sum_B = 0
    sum_F_uncertainty_squared = 0
    sum_B_uncertainty_squared = 0

    # Loop over the bins
    for bin_idx in range(1, num_bins + 1):
        bin_content = histogram.GetBinContent(bin_idx)
        bin_error = histogram.GetBinError(bin_idx)
        # Forward hemisphere (second half of bins)
        if bin_idx > half:  
            sum_F += bin_content
            sum_F_uncertainty_squared += bin_error**2
        # Backward hemisphere (first half of bins)
        else:  
            sum_B += bin_content
            sum_B_uncertainty_squared += bin_error**2

    # Compute FB asymmetry
    if sum_F + sum_B == 0:
        return 0, 0, 0, 0  # Avoid division by zero
    fb_asymmetry = (sum_F - sum_B) / (sum_F + sum_B)

    # Compute uncertainties
    sum_F_uncertainty = sum_F_uncertainty_squared**0.5
    sum_B_uncertainty = sum_B_uncertainty_squared**0.5

    fb_asymmetry_uncertainty = (
        2 * ((sum_B * sum_F_uncertainty)**2 + (sum_F * sum_B_uncertainty)**2)**0.5
        / (sum_F + sum_B)**2
    )

    # Compute coefficient and its uncertainty
    coefficient = fb_asymmetry * multiplier
    coefficient_uncertainty = fb_asymmetry_uncertainty * multiplier

    return fb_asymmetry, fb_asymmetry_uncertainty, coefficient, coefficient_uncertainty

def main():
    # Parse input arguments
    parser = argparse.ArgumentParser(description="Compute S.C. coefficients from histograms.")
    parser.add_argument("input_file", help="Path to the input ROOT file.")
    args = parser.parse_args()

    # Open the ROOT file
    root_file = ROOT.TFile.Open(args.input_file, "READ")
    if not root_file or root_file.IsZombie():
        print("Error: Could not open the ROOT file.")
        return

    # Open a text file to save the results with input file name
    # with open("extractedCoefficients.txt", "w") as output_file:
    with open(f"{args.input_file.split('/')[-1].replace('.root', '')}_extractedCoefficients.txt", "w") as output_file:
        # Loop over histograms in the dictionary
        for hist_name, multiplier in coeff_multipliers.items():
            histogram = root_file.Get(hist_name)
            if not histogram:
                print("Error: Histogram '{}' not found.".format(hist_name))
                continue

            # Compute FB asymmetry and coefficient
            fb_asymmetry, fb_asymmetry_uncertainty, coefficient, coefficient_uncertainty = compute_fb_asymmetry(histogram, multiplier)
            print("Histogram: {}".format(hist_name))
            print("  FB Asymmetry: {:.4f} +/- {:.4f}".format(fb_asymmetry, fb_asymmetry_uncertainty))
            print("  Coefficient: {:.4f} +/- {:.4f}".format(coefficient, coefficient_uncertainty))

            # Write the results to the text file
            output_file.write("Histogram: {}\n".format(hist_name))
            output_file.write("  FB Asymmetry: {:.4f} +/- {:.4f}\n".format(fb_asymmetry, fb_asymmetry_uncertainty))
            output_file.write("  Coefficient: {:.4f} +/- {:.4f}\n".format(coefficient, coefficient_uncertainty))
            output_file.write("\n")

    # Close the ROOT file
    root_file.Close()

if __name__ == "__main__":
    main()