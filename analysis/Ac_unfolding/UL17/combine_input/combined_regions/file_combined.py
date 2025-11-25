import ROOT
from optparse import OptionParser
import os
import sys

parser = OptionParser()
parser.add_option("-y", "--year", dest="year", help="Specify the year (UL18, UL17, preUL16, postUL16)", type='str')
parser.add_option("-m", "--mass_range", dest="mass_range", help="Specify the mass range (0_500, 500_750, 750-1000, 1000-1500, 1500Inf)", type='str')
parser.add_option("-l", "--lepton_flavor", dest="lepton_flavor", help="Specify the lepton flavor (electron, muon)", type='str')
parser.add_option("-r", "--region", dest="region", help="Specify the region (SR, CR1, CR2)", type='str')
parser.add_option("-s", "--region_score", dest="region_score", help="Specify the region score (0, 1, 2)", type='str')


(options, args) = parser.parse_args()

year = options.year if options.year else "UL17"  
mass_range = options.mass_range if options.mass_range else "0_500"  
lepton_flavor = options.lepton_flavor if options.lepton_flavor else "muon"


def copy_histograms(input_files, output_file_path):
    output_file = ROOT.TFile(output_file_path, 'RECREATE')

    directories = {
        '../individual_files/dY_{}_{}_{}_CR1.root'.format(year, lepton_flavor, mass_range) : 'CR1',
        '../individual_files/dY_{}_{}_{}_CR2.root'.format(year, lepton_flavor, mass_range): 'CR2',
        '../individual_files/dY_{}_{}_{}_SR.root'.format(year, lepton_flavor, mass_range) : 'SR'
    }

    for input_file_name, dir_name in directories.items():
        input_file = ROOT.TFile(input_file_name, 'READ')
        if input_file.IsOpen():
            print("Processing file: {}".format(input_file_name))

            output_dir = output_file.mkdir(dir_name)
            output_dir.cd()

            for key in input_file.GetListOfKeys():
                obj = key.ReadObj()
                if obj.InheritsFrom('TH1'):
                    hist = obj.Clone()
                    hist.Write()

            input_file.Close()
        else:
            print("Failed to open file: {}".format(input_file_name))

    output_file.Write()
    output_file.Close()

input_files = {
    '../individual_files/dY_{}_{}_{}_CR1.root'.format(year, lepton_flavor, mass_range) : 'CR1',
    '../individual_files/dY_{}_{}_{}_CR2.root'.format(year, lepton_flavor, mass_range): 'CR2',
    '../individual_files/dY_{}_{}_{}_SR.root'.format(year, lepton_flavor, mass_range) : 'SR'
}
output_file_path = 'dY_{}_{}_{}.root'.format(year, lepton_flavor, mass_range)

copy_histograms(input_files, output_file_path)
