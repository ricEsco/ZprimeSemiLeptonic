void CopyDirectory(TDirectory *source, TDirectory *target);

void salvage(const char* infile = "uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_UL17_23.root",
             const char* outfile = "salvaged3.root")
{
    // Open the corrupted file (READ only). If it's truly bad, even this might crash.
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open input file or file is zombie: " << infile << std::endl;
        return;
    }

    // Create a brand-new output file
    TFile *fout = new TFile(outfile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Cannot create output file: " << outfile << std::endl;
        fin->Close();
        return;
    }

    // Print statistics about the original file
    std::cout << "Original file size: " << fin->GetSize() << " bytes" << std::endl;

    // Loop over keys in the top-level directory of the input file
    TIter nextkey(fin->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        try {
            // Get class name
            TClass *cl = gROOT->GetClass(key->GetClassName());
            if (!cl) continue;

            // 1) If it's a directory (like TDirectoryFile), copy it *recursively*
            if (cl->InheritsFrom("TDirectory")) {
                TDirectory *sourceDir = (TDirectory*)key->ReadObj();
                if (!sourceDir) {
                    std::cerr << "Failed to read directory: " << key->GetName() << std::endl;
                    continue;
                }
                
                // Create a directory in the output file with the same name
                fout->cd();
                
                // Check if directory already exists
                TDirectory *targetDir = nullptr;
                if (fout->Get(sourceDir->GetName())) {
                    std::cout << "Directory " << sourceDir->GetName() << " already exists, using existing one" << std::endl;
                    targetDir = (TDirectory*)fout->Get(sourceDir->GetName());
                } else {
                    targetDir = fout->mkdir(sourceDir->GetName(), sourceDir->GetTitle());
                }
                
                if (targetDir) {
                    // Copy the contents
                    CopyDirectory(sourceDir, targetDir);
                } else {
                    std::cerr << "Failed to create directory: " << sourceDir->GetName() << std::endl;
                }
                delete sourceDir;
            }
            // 2) If it's a histogram (TH1), just read and write it
            else if (cl->InheritsFrom("TH1")) {
                TH1 *h = (TH1*)key->ReadObj();
                if (!h) {
                    std::cerr << "Failed to read histogram: " << key->GetName() << std::endl;
                    continue;
                }
                fout->cd();
                h->Write(); // saves histogram in the new file
                delete h;   // clean up memory
            }
            // 3) If it's a TTree, try to recover it
            else if (cl->InheritsFrom("TTree")) {
                std::cout << "Attempting to salvage TTree: " << key->GetName() << std::endl;
                
                try {
                    TTree *tree = (TTree*)key->ReadObj();
                    if (!tree) {
                        std::cerr << "Failed to read TTree: " << key->GetName() << std::endl;
                        continue;
                    }
                    
                    // Write the tree to the output file
                    fout->cd();
                    TTree *newtree = tree->CloneTree();
                    if (newtree) {
                        newtree->Write();
                        std::cout << "Successfully salvaged TTree: " << key->GetName() 
                                  << " with " << newtree->GetEntries() << " entries" << std::endl;
                    } else {
                        std::cerr << "Failed to clone TTree: " << key->GetName() << std::endl;
                    }
                    delete tree;
                } catch (std::exception &e) {
                    std::cerr << "Exception while salvaging TTree " << key->GetName() 
                              << ": " << e.what() << std::endl;
                } catch (...) {
                    std::cerr << "Unknown exception while salvaging TTree " << key->GetName() << std::endl;
                }
            }
            // 4) Try to salvage anything else by directly copying it
            else {
                std::cout << "Attempting to salvage object of class: " << key->GetClassName() 
                          << " named: " << key->GetName() << std::endl;
                
                try {
                    TObject *obj = key->ReadObj();
                    if (!obj) {
                        std::cerr << "Failed to read object: " << key->GetName() << std::endl;
                        continue;
                    }
                    
                    fout->cd();
                    obj->Write();
                    std::cout << "Successfully salvaged object: " << key->GetName() << std::endl;
                    delete obj;
                } catch (std::exception &e) {
                    std::cerr << "Exception while salvaging object " << key->GetName() 
                              << ": " << e.what() << std::endl;
                } catch (...) {
                    std::cerr << "Unknown exception while salvaging object " << key->GetName() << std::endl;
                }
            }
        } catch (std::exception &e) {
            std::cerr << "Exception caught while processing key " << key->GetName() 
                      << ": " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown exception caught while processing key " << key->GetName() << std::endl;
        }
    }

    fout->Close();
    fin->Close();
    
    // Get size of the salvaged file
    TFile *fcheck = TFile::Open(outfile, "READ");
    if (fcheck && !fcheck->IsZombie()) {
        std::cout << "Salvaged file size: " << fcheck->GetSize() << " bytes" << std::endl;
        fcheck->Close();
    }
    
    std::cout << "Salvage done. Output is in " << outfile << std::endl;
}

// A helper function to copy all keys from one directory to another, recursively.
void CopyDirectory(TDirectory *source, TDirectory *target)
{
    if (!source || !target) {
        std::cerr << "Null source or target directory in CopyDirectory" << std::endl;
        return;
    }

    // Save current directory so we can come back later
    TDirectory *savdir = gDirectory;

    // Go to the target directory
    target->cd();

    // Get list of keys in the source directory
    TIter nextkey(source->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        try {
            TClass *cl = gROOT->GetClass(key->GetClassName());
            if (!cl) continue;

            if (cl->InheritsFrom("TDirectory")) {
                // It's a subdirectory, so recurse
                TDirectory *subdir = (TDirectory*)key->ReadObj();
                if (!subdir) {
                    std::cerr << "Failed to read subdirectory: " << key->GetName() << std::endl;
                    continue;
                }
                
                // Check if directory already exists
                TDirectory *newSubDir = nullptr;
                if (target->Get(subdir->GetName())) {
                    std::cout << "Subdirectory " << subdir->GetName() << " already exists, using existing one" << std::endl;
                    newSubDir = (TDirectory*)target->Get(subdir->GetName());
                } else {
                    newSubDir = target->mkdir(subdir->GetName(), subdir->GetTitle());
                }
                
                if (newSubDir) {
                    CopyDirectory(subdir, newSubDir);
                } else {
                    std::cerr << "Failed to create subdirectory: " << subdir->GetName() << std::endl;
                }
                delete subdir;
            }
            else if (cl->InheritsFrom("TH1")) {
                // It's a histogram
                TH1 *h = (TH1*)key->ReadObj();
                if (!h) {
                    std::cerr << "Failed to read histogram: " << key->GetName() << std::endl;
                    continue;
                }
                h->Write();
                delete h;
            }
            // Also try to salvage TTrees within directories
            else if (cl->InheritsFrom("TTree")) {
                std::cout << "Attempting to salvage TTree in subdirectory: " << key->GetName() << std::endl;
                
                try {
                    TTree *tree = (TTree*)key->ReadObj();
                    if (!tree) {
                        std::cerr << "Failed to read TTree in subdirectory: " << key->GetName() << std::endl;
                        continue;
                    }
                    
                    // Write the tree to the target directory
                    TTree *newtree = tree->CloneTree();
                    if (newtree) {
                        newtree->Write();
                        std::cout << "Successfully salvaged TTree in subdirectory: " << key->GetName() 
                                  << " with " << newtree->GetEntries() << " entries" << std::endl;
                    } else {
                        std::cerr << "Failed to clone TTree in subdirectory: " << key->GetName() << std::endl;
                    }
                    delete tree;
                } catch (std::exception &e) {
                    std::cerr << "Exception while salvaging TTree " << key->GetName() 
                              << ": " << e.what() << std::endl;
                } catch (...) {
                    std::cerr << "Unknown exception while salvaging TTree " << key->GetName() << std::endl;
                }
            }
            else {
                // Try to salvage any other object type
                std::cout << "Attempting to salvage object in subdirectory, class: " << key->GetClassName() 
                          << " named: " << key->GetName() << std::endl;
                
                try {
                    TObject *obj = key->ReadObj();
                    if (!obj) {
                        std::cerr << "Failed to read object in subdirectory: " << key->GetName() << std::endl;
                        continue;
                    }
                    
                    obj->Write();
                    std::cout << "Successfully salvaged object in subdirectory: " << key->GetName() << std::endl;
                    delete obj;
                } catch (std::exception &e) {
                    std::cerr << "Exception while salvaging object " << key->GetName() 
                              << ": " << e.what() << std::endl;
                } catch (...) {
                    std::cerr << "Unknown exception while salvaging object " << key->GetName() << std::endl;
                }
            }
        } catch (std::exception &e) {
            std::cerr << "Exception caught while processing key " << key->GetName() 
                      << ": " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown exception caught while processing key " << key->GetName() << std::endl;
        }
    }

    // Go back to original directory
    if (savdir) savdir->cd();
}
