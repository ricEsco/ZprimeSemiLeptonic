#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TColor.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

void xtalk_test()
{
    gROOT->SetBatch(kTRUE);

    TFile* f_injtype1 = new TFile("Run000000_PixelAlive-13.root");
    TFile* f_injtype5 = new TFile("Run000001_PixelAlive-11.root");
    TFile* f_injtype6 = new TFile("Run000002_PixelAlive-15.root");

    // Efficiency thresholds
    float alive_eff     = 0.9;
    float coupled_eff   = 0.00001;
    float uncoupled_eff = 0.00001;

    // Base directory and object naming
    string baseDir      = "Detector/Board_0/OpticalGroup_0/Hybrid_0/Chip_";
    string shortBaseDir = "D_B(0)_O(0)_H(0)_";
    string module       = "SH0058"; 

    // --- Loop over all chips 12, 13, 14, 15 ----EDIT FOR DUAL
    for (int ch = 12; ch <= 15; ch++)
    {
        // Convert chip ID to string
        string chip = to_string(ch);

        TCanvas* c0_pixelalive1 = (TCanvas*) f_injtype1->Get(
            (baseDir + chip + "/" + shortBaseDir + "PixelAlive_Chip(" + chip + ")").c_str()
        );
        TH2F* h_pixelalive1 = c0_pixelalive1
            ? (TH2F*) c0_pixelalive1->GetPrimitive((shortBaseDir + "PixelAlive_Chip(" + chip + ")").c_str())
            : nullptr;

        TCanvas* c0_pixelalive5 = (TCanvas*) f_injtype5->Get(
            (baseDir + chip + "/" + shortBaseDir + "PixelAlive_Chip(" + chip + ")").c_str()
        );
        TH2F* h_pixelalive5 = c0_pixelalive5
            ? (TH2F*) c0_pixelalive5->GetPrimitive((shortBaseDir + "PixelAlive_Chip(" + chip + ")").c_str())
            : nullptr;

        TCanvas* c0_pixelalive6 = (TCanvas*) f_injtype6->Get(
            (baseDir + chip + "/" + shortBaseDir + "PixelAlive_Chip(" + chip + ")").c_str()
        );
        TH2F* h_pixelalive6 = c0_pixelalive6
            ? (TH2F*) c0_pixelalive6->GetPrimitive((shortBaseDir + "PixelAlive_Chip(" + chip + ")").c_str())
            : nullptr;

        // Check if any histogram is missing
        if (!h_pixelalive1 || !h_pixelalive5 || !h_pixelalive6) {
            cerr << "Error: Could not retrieve PixelAlive histograms for chip " << chip << endl;
            continue;
        }

        // Dimensions
        int nColumns = h_pixelalive1->GetXaxis()->GetNbins();
        int nRows    = h_pixelalive1->GetYaxis()->GetNbins();

        // Vectors to store pixel info
        vector<int> dead_row, dead_col;
        vector<int> confirmed_row, confirmed_col;

        // Clone the injection type 5 histogram to store missing bumps
        TH2F* h_confirmed2D = (TH2F*) h_pixelalive5->Clone(("h_confirmed2D_"+chip).c_str());
        h_confirmed2D->SetTitle(("Missing bump map for chip " + chip).c_str());

        // Analyze pixels
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nColumns; j++)
            {
                bool detectable = true;
                // 1) Check if it's alive in injection type 1, 90%
                if (h_pixelalive1->GetBinContent(j+1, i+1) < alive_eff) {
                    dead_row.push_back(i);
                    dead_col.push_back(j);
                    detectable = false;
                }

                // check injection types 5 & 6 thresholds
                if (detectable &&
                    h_pixelalive1->GetBinContent(j+1, i+1) >= alive_eff &&
                    h_pixelalive5->GetBinContent(j+1, i+1) <= coupled_eff &&
                    h_pixelalive6->GetBinContent(j+1, i+1) <= uncoupled_eff)
                {
                    // Missing bump
                    confirmed_row.push_back(i);
                    confirmed_col.push_back(j);
                    h_confirmed2D->SetBinContent(j+1, i+1, 1);
                }
                else {
                    // Good pixel
                    h_confirmed2D->SetBinContent(j+1, i+1, 0);
                }
            }
        }

        // Print summary
        cout << "Chip " << ch << ":\n"
             << "  alive_eff: " << alive_eff
             << ", coupled_eff: " << coupled_eff
             << ", uncoupled_eff: " << uncoupled_eff << endl;
        cout << "  dead pixels: " << dead_row.size() << endl;
        cout << "  missing bumps: " << confirmed_row.size() << endl;
        cout << "  Missing bump locations (row, col):\n";
        for (size_t k = 0; k < confirmed_row.size(); ++k) {
            cout << "    (" << confirmed_row[k] << ", " << confirmed_col[k] << ")\n";
        }

        // Save missing pixels to a file named with the chip ID
        ofstream outfile(("missing_crosstalk_pixels_" + chip + ".txt").c_str());
        for (size_t k = 0; k < confirmed_row.size(); ++k) {
            outfile << confirmed_row[k] << "," << confirmed_col[k] << endl;
        }
        outfile.close();

        // ---- 2-Color palette: 0 => yellow, 1 => black
        const Int_t nColors = 2;
        Int_t palette[nColors];
        palette[0] = TColor::GetColor("#FFFF00"); // Yellow for good pixels
        palette[1] = TColor::GetColor("#000000"); // Black for missing bumps

        gStyle->SetPalette(nColors, palette);
        h_confirmed2D->SetContour(nColors);
        h_confirmed2D->SetMinimum(-0.5);
        h_confirmed2D->SetMaximum(1.5);

        // Plot
        TCanvas* c_confirmed2D = new TCanvas(("c_confirmed2D_" + chip).c_str(),
                                             ("Missing bump map for chip " + chip).c_str(),
                                             800, 768);
        gStyle->SetOptStat(0);
        h_confirmed2D->Draw("colz");
        c_confirmed2D->SaveAs(("crosstalk_chip" + chip + ".png").c_str());
    } // end loop over ch=12..15
}
