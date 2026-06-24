from __future__ import print_function
import ROOT
from ttbargen_py import TTbarGenPy

fn = "/data/dust/user/ricardo/output_uhh2_Entanglement_Reco/UL18/preDNNselection/both/mergedFiles/uhh2.AnalysisModuleRunner.MC.TTToSemiLeptonic_muon.root"
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
while n_tau < 50 and i < 30000:
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
        if shown < 8:
            print("evt %d depth %d: tau(idx %d) -> daughters %s" % (i-1, depth, p.index(), [d.pdgId() for d in ds]))
        leps = [d for d in ds if abs(d.pdgId()) in (11, 13)]
        taus = [d for d in ds if abs(d.pdgId()) == 15]
        if leps:   found = leps[0]; break
        elif taus: p = taus[0]; depth += 1
        else:      break
    if found is not None:
        n_found += 1
        if shown < 8:
            print("   -> final lepton pdgId=%d pt=%.1f" % (found.pdgId(), found.v4().pt())); shown += 1
    elif shown < 8:
        print("   -> NO charged-lepton daughter in collection"); shown += 1

print("\ntauhad checked: %d | tau->e/mu found in gen record: %d (%.0f%%)" %
      (n_tau, n_found, 100.0*n_found/max(n_tau,1)))