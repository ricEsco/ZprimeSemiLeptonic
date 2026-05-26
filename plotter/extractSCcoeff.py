# -*- coding: utf-8 -*-
import sys
from ROOT import *

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
    sum_F = 0.0
    sum_B = 0.0
    sum_F_uncertainty_squared = 0.0
    sum_B_uncertainty_squared = 0.0

    n_bins = histogram.GetNbinsX()
    for i in range(1, n_bins + 1):
        bin_center = histogram.GetXaxis().GetBinCenter(i)
        bin_content = histogram.GetBinContent(i)
        bin_error = histogram.GetBinError(i)
        if bin_center >= 0:
            sum_F += bin_content
            sum_F_uncertainty_squared += bin_error**2
        else:
            sum_B += bin_content
            sum_B_uncertainty_squared += bin_error**2

    # Compute FB asymmetry
    if sum_F + sum_B == 0:
        return 0, 0, 0, 0  # Avoid division by zero
    fb_asymmetry = (sum_F - sum_B) / float(sum_F + sum_B)

    # Compute uncertainties
    sum_F_uncertainty = sum_F_uncertainty_squared**0.5
    sum_B_uncertainty = sum_B_uncertainty_squared**0.5

    fb_asymmetry_uncertainty = (
        2 * ((sum_B * sum_F_uncertainty)**2 + (sum_F * sum_B_uncertainty)**2)**0.5
        / float((sum_F + sum_B)**2)
    )

    # Compute coefficient and its uncertainty
    coefficient = fb_asymmetry * multiplier
    coefficient_uncertainty = fb_asymmetry_uncertainty * multiplier

    return fb_asymmetry, fb_asymmetry_uncertainty, coefficient, coefficient_uncertainty

def main():
    # Parse input arguments
    if len(sys.argv) < 2:
        print "Usage: python extractCoeff.py <input_file>"
        return

    input_file = sys.argv[1]

    root_file = ROOT.TFile.Open(input_file, "READ")
    if not root_file or root_file.IsZombie():
        print "Error: Could not open the ROOT file."
        return

    output_filename = "%s_extractedCoefficients.txt" % input_file.split('/')[-1].replace('.root', '')
    output_file = open(output_filename, "w")
    try:
        for hist_name, multiplier in coeff_multipliers.iteritems():
            histogram = root_file.Get(hist_name)
            if not histogram:
                print "Error: Histogram '%s' not found." % hist_name
                continue

            fb_asymmetry, fb_asymmetry_uncertainty, coefficient, coefficient_uncertainty = compute_fb_asymmetry(histogram, multiplier)
            print "Histogram: %s" % hist_name
            print "  FB Asymmetry: %.4f +/- %.4f" % (fb_asymmetry, fb_asymmetry_uncertainty)
            print "  Coefficient: %.4f +/- %.4f" % (coefficient, coefficient_uncertainty)

            output_file.write("Histogram: %s\n" % hist_name)
            output_file.write("  FB Asymmetry: %.4f +/- %.4f\n" % (fb_asymmetry, fb_asymmetry_uncertainty))
            output_file.write("  Coefficient: %.4f +/- %.4f\n" % (coefficient, coefficient_uncertainty))
            output_file.write("\n")
    finally:
        output_file.close()

    root_file.Close()

if __name__ == "__main__":
    main()