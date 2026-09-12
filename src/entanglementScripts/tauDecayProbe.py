from __future__ import print_function
import ROOT
from ttbargen_py import TTbarGenPy

fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_electron.root"
f = ROOT.TFile.Open(fn); t = f.Get("AnalysisTree")

def daughters_of(p, gps):
    out = []
    for j in range(gps.size()):
        gp = gps[j]
        m1 = gp.mother(gps, 1); m2 = gp.mother(gps, 2)
        if (m1 and m1.index() == p.index()) or (m2 and m2.index() == p.index()):
            out.append(gp)
    return out

n_tau = 0; n_found = 0; shown = 0; i = 0
while n_tau < 50000 and i < 100000:
    t.GetEntry(i); i += 1
    ttg = TTbarGenPy(t.GenParticles, False)
    if ttg.DecayChannel() != TTbarGenPy.e_tauhad:
        continue
    n_tau += 1
    gps = t.GenParticles
    p = ttg.ChargedLepton()           # the tau
    depth = 0; found = None
    while abs(p.pdgId()) == 15 and depth < 12:
        ds = daughters_of(p, gps)
        if shown < 20:
            print("evt %d depth %d: tau(idx %d) -> daughters %s" % (i-1, depth, p.index(), [d.pdgId() for d in ds]))
        leps = [d for d in ds if abs(d.pdgId()) in (11, 13)]
        taus = [d for d in ds if abs(d.pdgId()) == 15]
        if leps:   found = leps[0]; break
        elif taus: p = taus[0]; depth += 1
        else:      break
    if found is not None:
        n_found += 1
        if shown < 20:
            print("   -> final lepton pdgId=%d pt=%.1f" % (found.pdgId(), found.v4().pt())); shown += 1
    elif shown < 20:
        print("   -> NO charged-lepton daughter in collection"); shown += 1

print("\ntauhad checked: %d | tau->e/mu found in gen record: %d (%.0f%%)" %
      (n_tau, n_found, 100.0*n_found/max(n_tau,1)))

# Ran over 100000 events in the muon channel analysis output file, here is the output:
'''
evt 33 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 43 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 49 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 74 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 102 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 106 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 120 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 138 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 150 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 190 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 206 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 209 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 214 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 221 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 228 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 233 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 234 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 247 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 255 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 256 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection

tauhad checked: 7850 | tau->e/mu found in gen record: 0 (0%)
'''
# This shows us the tau decays are not present in the GenParticles collection

# same results for electron analysis output file:
'''
evt 0 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 21 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 34 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 42 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 58 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 73 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 79 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 91 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 92 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 114 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 121 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 155 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 160 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 176 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 178 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 179 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 220 depth 0: tau(idx 9) -> daughters []
   -> NO charged-lepton daughter in collection
evt 248 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 252 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection
evt 259 depth 0: tau(idx 11) -> daughters []
   -> NO charged-lepton daughter in collection

tauhad checked: 6415 | tau->e/mu found in gen record: 0 (0%)
'''