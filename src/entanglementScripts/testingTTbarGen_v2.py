from __future__ import print_function
import ROOT, glob, os

# 1) try to expose the C++ TTbarGen from compiled UHH2 libraries
libdir = os.path.join(os.environ["CMSSW_BASE"], "lib", os.environ["SCRAM_ARCH"])
cands = sorted(set(glob.glob(os.path.join(libdir, "*UHH2*ommon*")) +
                   glob.glob(os.path.join(libdir, "*UHH2*.so"))))
print("candidate UHH2 libs:")
for c in cands: print("   ", os.path.basename(c))
for c in cands:
    ROOT.gSystem.Load(c)                       # load errors here are OK; we check below
print("TTbarGen available after loads:", hasattr(ROOT, "TTbarGen"))

# open a file and grab one event's gen particles
FN = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root"
f = ROOT.TFile.Open(FN); t = f.Get("AnalysisTree"); t.GetEntry(0)
gps = t.GenParticles

if hasattr(ROOT, "TTbarGen"):
    ttg = ROOT.TTbarGen(gps, False)
    print("OK direct use -> decayChannel:", ttg.DecayChannel(), " isSL:", ttg.IsSemiLeptonicDecay())
    print("  BHad pt:", ttg.BHad().v4().pt() if ttg.IsSemiLeptonicDecay() else "n/a")
else:
    # 2) probe the GenParticle API so we can port TTbarGen faithfully
    p0 = gps[0]
    print("\nGenParticle methods:", [m for m in dir(p0) if not m.startswith("_")])
    for p in gps:
        if abs(p.pdgId()) == 6:
            print("\nprobing top pdgId", p.pdgId())
            for call in ["index", "daughter1", "daughter2", "mother1", "mother2", "status", "charge"]:
                try:    print("  %s() =" % call, getattr(p, call)())
                except Exception as e: print("  %s() ERR:" % call, e)
            try:
                d1 = p.daughter(gps, 1); d2 = p.daughter(gps, 2)
                print("  daughter(gps,1).pdgId =", d1.pdgId() if d1 else None)
                print("  daughter(gps,2).pdgId =", d2.pdgId() if d2 else None)
            except Exception as e:
                print("  daughter(gps,i) ERR:", e)
            break