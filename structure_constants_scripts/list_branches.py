# python list_branches.py EFT_700_900.root

import ROOT
import sys

def list_branches(filename, treename="AnalysisTree"):
    # List all branches in a ROOT file
    try:
        f = ROOT.TFile.Open(filename)
        if not f or f.IsZombie():
            print("Error: Could not open ROOT file: %s" % filename)
            return
        
        t = f.Get(treename)
        if not t:
            print("Error: Could not find TTree: %s" % treename)
            print("Available TTrees:")
            for key in f.GetListOfKeys():
                if key.GetClassName().startswith("TTree"):
                    print("  - %s" % key.GetName())
            return
        
        print("Branches in %s:%s:" % (filename, treename))
        for branch in t.GetListOfBranches():
            print("  - %s" % branch.GetName())
        
        f.Close()
    except Exception as e:
        print("Error: %s" % str(e))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: %s <root_file> [tree_name]" % sys.argv[0])
        sys.exit(1)
    
    filename = sys.argv[1]
    treename = sys.argv[2] if len(sys.argv) > 2 else "AnalysisTree"
    
    list_branches(filename, treename) 