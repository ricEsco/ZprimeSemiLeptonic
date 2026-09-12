from __future__ import print_function
import ROOT, glob, os

'''
This script is meant to investigate the use of the TTbarGen class from UHH2 independently of the rest of the UHH2 framework. 
TTbarGen is a C++ class that provides an interface to select the ttbar system particles at gen-level from Monte Carlo simulations:

uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/common/src/TTbarGen.cxx 

The script attempts to load the necessary UHH2 libraries and check if the TTbarGen class is available. 
If it is, it demonstrates how to use it to analyze a sample event. 
If not, it probes the GenParticle API to understand how to access the relevant information for top quark decays.
'''

### 1) Direct Access: try to expose the C++ TTbarGen from compiled UHH2 libraries ###
libdir = os.path.join(os.environ["CMSSW_BASE"], "lib", os.environ["SCRAM_ARCH"])
cands = sorted(set(glob.glob(os.path.join(libdir, "*UHH2*ommon*")) +
                   glob.glob(os.path.join(libdir, "*UHH2*.so"))))
print("candidate UHH2 libs:")
for c in cands: print("   ", os.path.basename(c))
for c in cands:
    ROOT.gSystem.Load(c)                       # load errors here are OK; we check below
print("TTbarGen available after loads:", hasattr(ROOT, "TTbarGen"))

# open a file and grab one event's gen particles
inFile = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root"
rootFile = ROOT.TFile.Open(inFile); rootTree = rootFile.Get("AnalysisTree"); rootTree.GetEntry(0)
genParticles = rootTree.GenParticles

if hasattr(ROOT, "TTbarGen"):
    ttg = ROOT.TTbarGen(genParticles, False)
    print("OK direct use -> decayChannel:", ttg.DecayChannel(), " isSL:", ttg.IsSemiLeptonicDecay())
    print("  BHad pt:", ttg.BHad().v4().pt() if ttg.IsSemiLeptonicDecay() else "n/a")

### 2) Indirect Access: probe the GenParticle API so we can port TTbarGen faithfully ###
else:
    particle0 = genParticles[0]
    print("\nGenParticle methods:", [m for m in dir(particle0) if not m.startswith("_")])
    for particle in genParticles:
        # select a top quark to probe its methods
        if abs(particle.pdgId()) == 6:
            print("\nprobing top pdgId", particle.pdgId())

            # probe the methods of the GenParticle class
            for call in ["index", "daughter1", "daughter2", "mother1", "mother2", "status", "charge"]:
                try:    print("  %s() =" % call, getattr(particle, call)())
                except Exception as e: print("  %s() ERR:" % call, e)

            # probe the daughter(genParticles, i) method which is used in TTbarGen
            try:
                d1 = particle.daughter(genParticles, 1); d2 = particle.daughter(genParticles, 2)
                print("  daughter(genParticles,1).pdgId =", d1.pdgId() if d1 else None)
                print("  daughter(genParticles,2).pdgId =", d2.pdgId() if d2 else None)
            except Exception as e:
                print("  daughter(genParticles,i) ERR:", e)
            break

'''
Results: 
Direct Accesss did not work. 
Indirect Access did work so we can at least use the same GenParticle methods to replicate the module's logic in selecting the ttbar system particles.
Below is the probe output:

GenParticle methods: ['charge', 'daughter', 'daughter1', 'daughter2', 'energy', 'eta', 'hadronFlavour', 'index', 'mother', 'mother1', 'mother2', 'partonFlavour', 'pdgId', 'phi', 'pt', 'set_charge', 'set_daughter1', 'set_daughter2', 'set_energy', 'set_eta', 'set_hadronFlavour', 'set_index', 'set_mother1', 'set_mother2', 'set_partonFlavour', 'set_pdgId', 'set_phi', 'set_pt', 'set_spin', 'set_status', 'set_v4', 'spin', 'status', 'v4']

probing top pdgId 6
  index() = 2
  daughter1() = 5
  daughter2() = 6
  mother1() = 0
  mother2() = 1
  status() = 22
  charge() = 0
  daughter(genParticles,1).pdgId = 24
  daughter(genParticles,2).pdgId = 5
Error in <THashList::Delete>: A list is accessing an object (0xf0a3e30) already deleted (list name = THashList)
*** Error in `python': double free or corruption (!prev): 0x0000000005ececa0 ***

'''