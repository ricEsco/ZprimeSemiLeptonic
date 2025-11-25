import ROOT
import sys
import array

def crosstalk(moduleName="SH0054"):
    ROOT.gROOT.SetBatch(True)
    # Turn off all stats boxes
    # ROOT.gStyle.SetOptStat(0)

    # Determine chip list based on module name prefix
    if moduleName.startswith("SH"):
        chipList = [12, 13, 14, 15]  # Quad configuration
    elif moduleName.startswith("RH"):
        chipList = [12, 13]          # Dual configuration
    else:
        print(f"Unrecognized module name prefix: {moduleName}")
        return

    # Open your input files (Edit the filenames!)
    f_injtype1 = ROOT.TFile("Run000000_PixelAlive-13.root")  # injection type 1
    f_injtype5 = ROOT.TFile("Run000001_PixelAlive-11.root")  # injection type 5
    f_injtype6 = ROOT.TFile("Run000002_PixelAlive-15.root")  # injection type 6

    # Create an output ROOT file to store all results
    f_out = ROOT.TFile("crosstalk_analysis_output.root", "RECREATE")

    # Efficiency thresholds
    alive_eff     = 0.90
    coupled_eff   = 1e-5
    uncoupled_eff = 1e-5

    # Base directory and object naming
    baseDir      = "Detector/Board_0/OpticalGroup_0/Hybrid_0/Chip_"
    shortBaseDir = "D_B(0)_O(0)_H(0)_"

    # Loop over the chips from chipList
    for ch in chipList:
        chipStr = str(ch)

        # Retrieve PixelAlive canvases and histograms for each injection type
        c0_pixelalive1 = f_injtype1.Get(f"{baseDir}{chipStr}/{shortBaseDir}PixelAlive_Chip({chipStr})")
        h_pixelalive1 = (c0_pixelalive1.GetPrimitive(f"{shortBaseDir}PixelAlive_Chip({chipStr})")
                         if c0_pixelalive1 else None)

        c0_pixelalive5 = f_injtype5.Get(f"{baseDir}{chipStr}/{shortBaseDir}PixelAlive_Chip({chipStr})")
        h_pixelalive5 = (c0_pixelalive5.GetPrimitive(f"{shortBaseDir}PixelAlive_Chip({chipStr})")
                         if c0_pixelalive5 else None)

        c0_pixelalive6 = f_injtype6.Get(f"{baseDir}{chipStr}/{shortBaseDir}PixelAlive_Chip({chipStr})")
        h_pixelalive6 = (c0_pixelalive6.GetPrimitive(f"{shortBaseDir}PixelAlive_Chip({chipStr})")
                         if c0_pixelalive6 else None)

        # Check if any histogram is missing
        if not (h_pixelalive1 and h_pixelalive5 and h_pixelalive6):
            print(f"Error: Could not retrieve PixelAlive histograms for chip {chipStr}")
            continue

        # Dimensions of the histogram
        nColumns = h_pixelalive1.GetXaxis().GetNbins()
        nRows    = h_pixelalive1.GetYaxis().GetNbins()

        # Vectors to store pixel information
        confirmed_row = []
        confirmed_col = []

        # Clone one histogram to serve as a 2D map of open bumps
        h_confirmed2D = h_pixelalive5.Clone(f"h_confirmed2D_{chipStr}")
        h_confirmed2D.SetTitle(f"Open bump map for chip {chipStr}")
        # Remove the stat box from this histogram
        h_confirmed2D.SetStats(False)

        # Analyze pixels: loop over all rows (i) and columns (j)
        for i in range(nRows):
            for j in range(nColumns):
                # 1) Check if the pixel is "alive" in injection type 1 (≥ 90% efficiency)
                alive = (h_pixelalive1.GetBinContent(j+1, i+1) >= alive_eff)

                # 2) For pixels that are alive, check if they show no crosstalk in injtype5 & injtype6
                noCrosstalk = (h_pixelalive5.GetBinContent(j+1, i+1) <= coupled_eff and
                               h_pixelalive6.GetBinContent(j+1, i+1) <= uncoupled_eff)

                if alive and noCrosstalk:
                    confirmed_row.append(i)
                    confirmed_col.append(j)
                    h_confirmed2D.SetBinContent(j+1, i+1, 1)
                else:
                    h_confirmed2D.SetBinContent(j+1, i+1, 0)

        # Print a summary for each chip
        print(f"Chip {chipStr}:")
        print(f"  Found {len(confirmed_row)} open bumps.")
        for idx in range(len(confirmed_row)):
            print(f"    ({confirmed_row[idx]}, {confirmed_col[idx]})")

        # Create a TTree for the numeric data
        t_openbump = ROOT.TTree(f"openbumpCoordinates_Chip{chipStr}", "Open bump coordinates")
        openbump_row = array.array('i', [0])
        openbump_col = array.array('i', [0])
        t_openbump.Branch("row", openbump_row, "row/I")
        t_openbump.Branch("col", openbump_col, "col/I")

        for idx in range(len(confirmed_row)):
            openbump_row[0] = confirmed_row[idx]
            openbump_col[0] = confirmed_col[idx]
            t_openbump.Fill()

        # Modify the 2D map style: 2-color palette
        nColors = 2
        palette = [0] * nColors
        palette[0] = ROOT.TColor.GetColor("#FFFF00")  # Yellow
        palette[1] = ROOT.TColor.GetColor("#000000")  # Black

        # Convert to array for SetPalette
        palette_arr = array.array('i', palette)
        ROOT.gStyle.SetPalette(nColors, palette_arr)
        h_confirmed2D.SetContour(nColors)
        h_confirmed2D.SetMinimum(-0.5)
        h_confirmed2D.SetMaximum(1.5)

        f_out.cd()
        h_confirmed2D.Write()  # 2D map
        t_openbump.Write()      # TTree with row & col

    f_out.Close()


# Optional: Take the module name from command-line arguments
if __name__ == "__main__":
    moduleName = "SH0054"
    if len(sys.argv) > 1:
        moduleName = sys.argv[1]
    crosstalk(moduleName)
