# ttbargen_py.py  --  pyROOT port of UHH2 TTbarGen.cxx
# Build per event from t.GenParticles; use it before the next GetEntry (proxies point into the entry).
from __future__ import print_function

class TTbarGenPy(object):
    # decay-channel enum (same order/meaning as TTbarGen.h)
    (e_had, e_ehad, e_muhad, e_tauhad, e_ee, e_mumu,
     e_tautau, e_emu, e_etau, e_mutau, e_notfound) = range(11)

    def __init__(self, genparticles, throw_on_failure=False):
        self.gps = genparticles
        self.m_type = self.e_notfound
        self.m_Top = self.m_Antitop = None
        self.m_WTop = self.m_WAntitop = None
        self.m_bTop = self.m_bAntitop = None
        self.m_Wdecay1 = self.m_Wdecay2 = None
        self.m_WMinusdecay1 = self.m_WMinusdecay2 = None
        self._ok = self._build(bool(throw_on_failure))

    def _has_top_mother(self, gp, top):
        m1 = gp.mother(self.gps, 1)
        m2 = gp.mother(self.gps, 2)
        return (bool(m1) and m1.index() == top.index()) or \
               (bool(m2) and m2.index() == top.index())

    def _build(self, throw):
        gps = self.gps
        n_top = n_antitop = 0
        top = antitop = None
        for i in range(gps.size()):
            genp = gps[i]
            if abs(genp.pdgId()) != 6:
                continue
            if genp.pdgId() == 6  and top is None:     top = genp
            if genp.pdgId() == -6 and antitop is None: antitop = genp

            w = genp.daughter(gps, 1)
            b = genp.daughter(gps, 2)
            if (not w) or (not b):
                if throw: raise RuntimeError("top has not ==2 daughters")
                return False
            if abs(w.pdgId()) != 24:                 # ensure w is the W, b the quark
                w, b = b, w
            if abs(w.pdgId()) != 24:                 # radiation workaround: find W by top-mother
                for j in range(gps.size()):
                    gp = gps[j]
                    if abs(gp.pdgId()) == 24 and self._has_top_mother(gp, genp):
                        w = gp; break
            if abs(b.pdgId()) not in (5, 3, 1):      # same workaround for the b
                for j in range(gps.size()):
                    gp = gps[j]
                    if abs(gp.pdgId()) in (5, 3, 1) and self._has_top_mother(gp, genp):
                        b = gp; break
            if abs(b.pdgId()) not in (5, 3, 1):
                if throw: raise RuntimeError("top has no b daughter")
                return False

            wd1 = w.daughter(gps, 1)
            wd2 = w.daughter(gps, 2)
            n_wdau = 0
            while n_wdau != 2:                       # follow W -> W copies until 2 daughters
                if wd1 and not wd2:
                    w = wd1
                    wd1 = w.daughter(gps, 1); wd2 = w.daughter(gps, 2)
                elif wd1 and wd2:
                    n_wdau = 2
                else:
                    if throw: raise RuntimeError("W has no daughters")
                    return False

            if genp.pdgId() == 6:
                self.m_Top, self.m_WTop, self.m_bTop = top, w, b
                self.m_Wdecay1, self.m_Wdecay2 = wd1, wd2
                n_top += 1
            else:
                self.m_Antitop, self.m_WAntitop, self.m_bAntitop = antitop, w, b
                self.m_WMinusdecay1, self.m_WMinusdecay2 = wd1, wd2
                n_antitop += 1

        if n_top != 1 or n_antitop != 1:
            if throw: raise RuntimeError("did not find exactly one top and one antitop")
            return False

        n_e = n_m = n_t = 0
        for wd in (self.m_Wdecay1, self.m_Wdecay2, self.m_WMinusdecay1, self.m_WMinusdecay2):
            a = abs(wd.pdgId())
            if   a == 11: n_e += 1
            elif a == 13: n_m += 1
            elif a == 15: n_t += 1
        if   n_e == 2:              self.m_type = self.e_ee
        elif n_e == 1 and n_m == 1: self.m_type = self.e_emu
        elif n_e == 1 and n_t == 1: self.m_type = self.e_etau
        elif n_m == 2:              self.m_type = self.e_mumu
        elif n_m == 1 and n_t == 1: self.m_type = self.e_mutau
        elif n_t == 2:              self.m_type = self.e_tautau
        elif n_e == 1:              self.m_type = self.e_ehad
        elif n_m == 1:              self.m_type = self.e_muhad
        elif n_t == 1:              self.m_type = self.e_tauhad
        else:                       self.m_type = self.e_had
        return True

    # ---- accessors (mirror TTbarGen.h) ----
    def DecayChannel(self):        return self.m_type
    def IsSemiLeptonicDecay(self): return self.m_type in (self.e_ehad, self.e_muhad, self.e_tauhad)
    def Top(self):      return self.m_Top
    def Antitop(self):  return self.m_Antitop
    def bTop(self):     return self.m_bTop
    def bAntitop(self): return self.m_bAntitop

    def _lep(self, gp): return abs(gp.pdgId()) in (11, 13, 15)
    def _nu(self, gp):  return abs(gp.pdgId()) in (12, 14, 16)
    def ChargedLepton(self):
        for wd in (self.m_Wdecay1, self.m_Wdecay2, self.m_WMinusdecay1, self.m_WMinusdecay2):
            if self._lep(wd): return wd
        raise RuntimeError("ChargedLepton: not an l+jets event")
    def Neutrino(self):
        for wd in (self.m_Wdecay1, self.m_Wdecay2, self.m_WMinusdecay1, self.m_WMinusdecay2):
            if self._nu(wd): return wd
        raise RuntimeError("Neutrino: not an l+jets event")
    def TopLep(self): return self.Top() if self.ChargedLepton().charge() > 0 else self.Antitop()
    def TopHad(self): return self.Top() if self.ChargedLepton().charge() < 0 else self.Antitop()
    def BLep(self):   return self.bTop() if self.ChargedLepton().charge() > 0 else self.bAntitop()
    def BHad(self):   return self.bTop() if self.ChargedLepton().charge() < 0 else self.bAntitop()