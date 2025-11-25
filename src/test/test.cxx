// --- TEMPLATE METHOD FILLING (systematics) ---
if(is_mc && is_tt){
    // 1) reco xi from the event-level DeltaY you already compute
    const float xi_reco = std::tanh(dyreco);
  
    // 2) need xi_gen present for lookup
    if(use_noac_evtweights_ && event.is_valid(h_xi_gen)){
      const double xi_gen_evt = event.get(h_xi_gen);
  
      // 3) fill the nominal-per-f templates
      for(const float fv : f_values){
        auto itW = noac_weights_map.find(fv);
        if(itW != noac_weights_map.end() && itW->second){
          const double w_noac = lookup_noac_weight(itW->second.get(), xi_gen_evt);
          std::string f_str = std::to_string(fv);
          std::replace(f_str.begin(), f_str.end(), '.', 'p');
          std::replace(f_str.begin(), f_str.end(), '-', 'm');
          const std::string name = std::string("DeltaY_xi_reco_6_f_") + f_str;
  
          auto itH = h_deltaY_xi_reco_map.find(name);
          if(itH != h_deltaY_xi_reco_map.end() && itH->second){
            itH->second->Fill(xi_reco, weight * w_noac);
          }
        }
      }
  
      // 4) optional: NoAC nuisance on nominal xi (choose your step, e.g. +/−0.2 around f=1)
      // If you want to follow your doc exactly, replace 0.2 with that prescription.
      const double f_step = 0.2;
      auto makeW = [&](double f)->const TH1D*{
        // reuse prepared ones when possible, else build on the fly from the summed gen
        auto it = noac_weights_map.find(f);
        if(it != noac_weights_map.end()) return it->second.get();
        return nullptr;
      };
      if(auto W_up = makeW(1.0 + f_step)){
        const double w_up = lookup_noac_weight(W_up, xi_gen_evt);
        DeltaY_xi_reco_6_NoAC_up->Fill(xi_reco, weight * w_up);
      }
      if(auto W_dn = makeW(1.0 - f_step)){
        const double w_dn = lookup_noac_weight(W_dn, xi_gen_evt);
        DeltaY_xi_reco_6_NoAC_down->Fill(xi_reco, weight * w_dn);
      }
    }
  
    // Also fill the nominal (no template) histogram
    DeltaY_xi_reco_6->Fill(xi_reco, weight);
  }
  