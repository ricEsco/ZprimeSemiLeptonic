# to test retrieving structure constants from ROOT file and calculating weights as well as comparing to the weights in the root file

import ROOT
import numpy as np
import sys
import argparse

def get_wc_names():
    return {
    0: "ctGRe",
    1: "ctGIm",
    2: "cQj18",
    3: "cQj38",
    4: "cQj11",
    5: "cQj31",
    6: "ctu8",
    7: "ctd8",
    8: "ctj8",
    9: "cQu8",
    10: "cQd8",
    11: "ctu1",
    12: "ctd1",
    13: "ctj1",
    14: "cQu1",
    15: "cQd1"
}

def get_wc_index(wc_name):
    # Get the index of a Wilson coefficient by name
    wc_names = get_wc_names()
    for idx, name in wc_names.items():
        if name == wc_name:
            return idx
    return -1  # Not found

def create_single_wc_scenario(wc_name, value=1.0):
    # a scenario with a single WC set to a value
    wc_values = [0.0] * 16
    idx = get_wc_index(wc_name)
    if idx >= 0:
        wc_values[idx] = value
    return wc_values

# Dictionary mapping pairs of WC indices to their weight indices
KNOWN_PAIRS = {
    # ctGRe (0) pairs: 15 pairs
    (0,1): 220,   # ctGRe_ctGIm_1_1 -> EFTrwgt218
    (0,2): 221,   # ctGRe_cQj18_1_1 -> EFTrwgt219
    (0,3): 222,   # ctGRe_cQj38_1_1 -> EFTrwgt220
    (0,4): 223,   # ctGRe_cQj11_1_1 -> EFTrwgt221
    (0,5): 224,   # ctGRe_cQj31_1_1 -> EFTrwgt222
    (0,6): 225,   # ctGRe_ctu8_1_1 -> EFTrwgt223
    (0,7): 226,   # ctGRe_ctd8_1_1 -> EFTrwgt224
    (0,8): 227,   # ctGRe_ctj8_1_1 -> EFTrwgt225
    (0,9): 228,   # ctGRe_cQu8_1_1 -> EFTrwgt226
    (0,10): 229,  # ctGRe_cQd8_1_1 -> EFTrwgt227
    (0,11): 230,  # ctGRe_ctu1_1_1 -> EFTrwgt228
    (0,12): 231,  # ctGRe_ctd1_1_1 -> EFTrwgt229
    (0,13): 232,  # ctGRe_ctj1_1_1 -> EFTrwgt230
    (0,14): 233,  # ctGRe_cQu1_1_1 -> EFTrwgt231
    (0,15): 234,  # ctGRe_cQd1_1_1 -> EFTrwgt232

    # ctGIm (1) pairs: 14 pairs
    (1,2): 236,   # ctGIm_cQj18_1_1 -> EFTrwgt234
    (1,3): 237,   # ctGIm_cQj38_1_1 -> EFTrwgt235
    (1,4): 238,   # ctGIm_cQj11_1_1 -> EFTrwgt236
    (1,5): 239,   # ctGIm_cQj31_1_1 -> EFTrwgt237
    (1,6): 240,   # ctGIm_ctu8_1_1 -> EFTrwgt238
    (1,7): 241,   # ctGIm_ctd8_1_1 -> EFTrwgt239
    (1,8): 242,   # ctGIm_ctj8_1_1 -> EFTrwgt240
    (1,9): 243,   # ctGIm_cQu8_1_1 -> EFTrwgt241
    (1,10): 244,  # ctGIm_cQd8_1_1 -> EFTrwgt242
    (1,11): 245,  # ctGIm_ctu1_1_1 -> EFTrwgt243
    (1,12): 246,  # ctGIm_ctd1_1_1 -> EFTrwgt244
    (1,13): 247,  # ctGIm_ctj1_1_1 -> EFTrwgt245
    (1,14): 248,  # ctGIm_cQu1_1_1 -> EFTrwgt246
    (1,15): 249,  # ctGIm_cQd1_1_1 -> EFTrwgt247

    # cQj18 (2) pairs: 13 pairs
    (2,3): 251,   # cQj18_cQj38_1_1 -> EFTrwgt249
    (2,4): 252,   # cQj18_cQj11_1_1 -> EFTrwgt250
    (2,5): 253,   # cQj18_cQj31_1_1 -> EFTrwgt251
    (2,6): 254,   # cQj18_ctu8_1_1 -> EFTrwgt252
    (2,7): 255,   # cQj18_ctd8_1_1 -> EFTrwgt253
    (2,8): 256,   # cQj18_ctj8_1_1 -> EFTrwgt254
    (2,9): 257,   # cQj18_cQu8_1_1 -> EFTrwgt255
    (2,10): 258,  # cQj18_cQd8_1_1 -> EFTrwgt256
    (2,11): 259,  # cQj18_ctu1_1_1 -> EFTrwgt257
    (2,12): 260,  # cQj18_ctd1_1_1 -> EFTrwgt258
    (2,13): 261,  # cQj18_ctj1_1_1 -> EFTrwgt259
    (2,14): 262,  # cQj18_cQu1_1_1 -> EFTrwgt260
    (2,15): 263,  # cQj18_cQd1_1_1 -> EFTrwgt261

    # cQj38 (3) pairs: 12 pairs
    (3,4): 265,   # cQj38_cQj11_1_1 -> EFTrwgt263
    (3,5): 266,   # cQj38_cQj31_1_1 -> EFTrwgt264
    (3,6): 267,   # cQj38_ctu8_1_1 -> EFTrwgt265
    (3,7): 268,   # cQj38_ctd8_1_1 -> EFTrwgt266
    (3,8): 269,   # cQj38_ctj8_1_1 -> EFTrwgt267
    (3,9): 270,   # cQj38_cQu8_1_1 -> EFTrwgt268
    (3,10): 271,  # cQj38_cQd8_1_1 -> EFTrwgt269
    (3,11): 272,  # cQj38_ctu1_1_1 -> EFTrwgt270
    (3,12): 273,  # cQj38_ctd1_1_1 -> EFTrwgt271
    (3,13): 274,  # cQj38_ctj1_1_1 -> EFTrwgt272
    (3,14): 275,  # cQj38_cQu1_1_1 -> EFTrwgt273
    (3,15): 276,  # cQj38_cQd1_1_1 -> EFTrwgt274

    # cQj11 (4) pairs: 11 pairs
    (4,5): 278,   # cQj11_cQj31_1_1 -> EFTrwgt276
    (4,6): 279,   # cQj11_ctu8_1_1 -> EFTrwgt277
    (4,7): 280,   # cQj11_ctd8_1_1 -> EFTrwgt278
    (4,8): 281,   # cQj11_ctj8_1_1 -> EFTrwgt279
    (4,9): 282,   # cQj11_cQu8_1_1 -> EFTrwgt280
    (4,10): 283,  # cQj11_cQd8_1_1 -> EFTrwgt281
    (4,11): 284,  # cQj11_ctu1_1_1 -> EFTrwgt282
    (4,12): 285,  # cQj11_ctd1_1_1 -> EFTrwgt283
    (4,13): 286,  # cQj11_ctj1_1_1 -> EFTrwgt284
    (4,14): 287,  # cQj11_cQu1_1_1 -> EFTrwgt285
    (4,15): 288,  # cQj11_cQd1_1_1 -> EFTrwgt286

    # cQj31 (5) pairs: 10 pairs
    (5,6): 290,   # cQj31_ctu8_1_1 -> EFTrwgt288
    (5,7): 291,   # cQj31_ctd8_1_1 -> EFTrwgt289
    (5,8): 292,   # cQj31_ctj8_1_1 -> EFTrwgt290
    (5,9): 293,   # cQj31_cQu8_1_1 -> EFTrwgt291
    (5,10): 294,  # cQj31_cQd8_1_1 -> EFTrwgt292
    (5,11): 295,  # cQj31_ctu1_1_1 -> EFTrwgt293
    (5,12): 296,  # cQj31_ctd1_1_1 -> EFTrwgt294
    (5,13): 297,  # cQj31_ctj1_1_1 -> EFTrwgt295
    (5,14): 298,  # cQj31_cQu1_1_1 -> EFTrwgt296
    (5,15): 299,  # cQj31_cQd1_1_1 -> EFTrwgt297

    # ctu8 (6) pairs: 9 pairs
    (6,7): 301,   # ctu8_ctd8_1_1 -> EFTrwgt299
    (6,8): 302,   # ctu8_ctj8_1_1 -> EFTrwgt300
    (6,9): 303,   # ctu8_cQu8_1_1 -> EFTrwgt301
    (6,10): 304,  # ctu8_cQd8_1_1 -> EFTrwgt302
    (6,11): 305,  # ctu8_ctu1_1_1 -> EFTrwgt303
    (6,12): 306,  # ctu8_ctd1_1_1 -> EFTrwgt304
    (6,13): 307,  # ctu8_ctj1_1_1 -> EFTrwgt305
    (6,14): 308,  # ctu8_cQu1_1_1 -> EFTrwgt306
    (6,15): 309,  # ctu8_cQd1_1_1 -> EFTrwgt307

    # ctd8 (7) pairs: 8 pairs
    (7,8): 311,   # ctd8_ctj8_1_1 -> EFTrwgt309
    (7,9): 312,   # ctd8_cQu8_1_1 -> EFTrwgt310
    (7,10): 313,  # ctd8_cQd8_1_1 -> EFTrwgt311
    (7,11): 314,  # ctd8_ctu1_1_1 -> EFTrwgt312
    (7,12): 315,  # ctd8_ctd1_1_1 -> EFTrwgt313
    (7,13): 316,  # ctd8_ctj1_1_1 -> EFTrwgt314
    (7,14): 317,  # ctd8_cQu1_1_1 -> EFTrwgt315
    (7,15): 318,  # ctd8_cQd1_1_1 -> EFTrwgt316

    # ctj8 (8) pairs: 7 pairs
    (8,9): 320,   # ctj8_cQu8_1_1 -> EFTrwgt318
    (8,10): 321,  # ctj8_cQd8_1_1 -> EFTrwgt319
    (8,11): 322,  # ctj8_ctu1_1_1 -> EFTrwgt320
    (8,12): 323,  # ctj8_ctd1_1_1 -> EFTrwgt321
    (8,13): 324,  # ctj8_ctj1_1_1 -> EFTrwgt322
    (8,14): 325,  # ctj8_cQu1_1_1 -> EFTrwgt323
    (8,15): 326,  # ctj8_cQd1_1_1 -> EFTrwgt324

    # cQu8 (9) pairs: 6 pairs
    (9,10): 328,  # cQu8_cQd8_1_1 -> EFTrwgt326
    (9,11): 329,  # cQu8_ctu1_1_1 -> EFTrwgt327
    (9,12): 330,  # cQu8_ctd1_1_1 -> EFTrwgt328
    (9,13): 331,  # cQu8_ctj1_1_1 -> EFTrwgt329
    (9,14): 332,  # cQu8_cQu1_1_1 -> EFTrwgt330
    (9,15): 333,  # cQu8_cQd1_1_1 -> EFTrwgt331

    # cQd8 (10) pairs: 5 pairs
    (10,11): 335, # cQd8_ctu1_1_1 -> EFTrwgt333
    (10,12): 336, # cQd8_ctd1_1_1 -> EFTrwgt334
    (10,13): 337, # cQd8_ctj1_1_1 -> EFTrwgt335
    (10,14): 338, # cQd8_cQu1_1_1 -> EFTrwgt336
    (10,15): 339, # cQd8_cQd1_1_1 -> EFTrwgt337

    # ctu1 (11) pairs: 4 pairs
    (11,12): 341, # ctu1_ctd1_1_1 -> EFTrwgt339
    (11,13): 342, # ctu1_ctj1_1_1 -> EFTrwgt340
    (11,14): 343, # ctu1_cQu1_1_1 -> EFTrwgt341
    (11,15): 344, # ctu1_cQd1_1_1 -> EFTrwgt342

    # ctd1 (12) pairs: 3 pairs
    (12,13): 346, # ctd1_ctj1_1_1 -> EFTrwgt344
    (12,14): 347, # ctd1_cQu1_1_1 -> EFTrwgt345
    (12,15): 348, # ctd1_cQd1_1_1 -> EFTrwgt346

    # ctj1 (13) pairs: 2 pairs
    (13,14): 350, # ctj1_cQu1_1_1 -> EFTrwgt348
    (13,15): 351, # ctj1_cQd1_1_1 -> EFTrwgt349

    # cQu1 (14) pairs: 1 pair
    (14,15): 353  # cQu1_cQd1_1_1 -> EFTrwgt351
}

# this function gets the index of the weight in the systweights array for a given WC configuration
# this is for the saved weights in the root file in order to compare to the calculated weights
def get_weight_index(wc_values):
    # Args:
    # wc_values: list of Wilson coefficient values (length: num_WCs)
    # Returns:
    # weight_idx: the index in the systweights array, or -1 if not found
    
    # SM point (all zeros)
    if all(v == 0.0 for v in wc_values):
        return 201  # SM point
    
    # Count nonzero WCs
    non_zero_indices = [i for i, v in enumerate(wc_values) if v != 0.0]
    num_non_zero = len(non_zero_indices)
    
    # Single WC case
    if num_non_zero == 1:
        idx = non_zero_indices[0]
        value = wc_values[idx]
        
        # Single WC = 1
        if value == 1.0:
            return 203 + idx  # Single WC=1 indices start at 203
        
        # Single WC = 2
        if value == 2.0:
            # Mapping for single=2
            single2_mapping = {
                0: 219,   # ctGRe
                1: 235,   # ctGIm
                2: 250,   # cQj18
                3: 264,   # cQj38
                4: 277,   # cQj11
                5: 289,   # cQj31
                6: 300,   # ctu8
                7: 310,   # ctd8
                8: 319,   # ctj8
                9: 327,   # cQu8
                10: 335,  # cQd8
                11: 340,  # ctu1
                12: 345,  # ctd1
                13: 349,  # ctj1
                14: 352,  # cQu1
                15: 354   # cQd1
            }
            return single2_mapping.get(idx, -1)
    
    # Pair WC case (both = 1)
    if num_non_zero == 2:
        idx1, idx2 = non_zero_indices
        val1, val2 = wc_values[idx1], wc_values[idx2]
        
        # Both values must be 1.0
        if val1 == 1.0 and val2 == 1.0:
            # idx1 < idx2 for the mapping
            if idx1 > idx2:
                idx1, idx2 = idx2, idx1
            
            # Look up in KNOWN_PAIRS
            return KNOWN_PAIRS.get((idx1, idx2), -1)
    
    # Other configurations not supported
    return -1

# this function calculates the weight using the structure constants and the given WC values
def calculate_weights(structure_constants, wc_values):
    # structure_constants: numpy array of structure constants (shape: poly_dim)
    # wc_values: list of WC values (length: num_WCs)
    
    # Returns:
    # total_weight: the calculated weight
    num_WCs = len(wc_values)

    # general order is: index 0 is the constant term, then linear terms (next num_WCs terms), then quadratic terms (next num_WCs terms), then cross terms (next num_WCs*(num_WCs-1)/2 terms)
    
    # Constant term (SM)
    c0 = structure_constants[0]
    
    # Linear terms
    # linear_terms: numpy array of linear terms (shape: num_WCs), structure constants
    linear_terms = structure_constants[1:1+num_WCs]
    # w_linear: the calculated weight for the linear terms
    w_linear = np.dot(linear_terms, wc_values)
    

    # Quadratic terms
    # quad_list: list of quadratic terms (shape: num_WCs), structure constants
    quad_list = []
    # Diagonal terms (WC^2)
    for i in range(num_WCs):
        quad_list.append(wc_values[i]**2)
    
    # Cross terms (WC_i * WC_j)
    for i in range(num_WCs):
        for j in range(i+1, num_WCs):
            quad_list.append(wc_values[i] * wc_values[j])
    
    # quad_terms: numpy array of quadratic terms (shape: num_WCs), structure constants
    idx_quad_start = 1 + num_WCs
    quad_terms = structure_constants[idx_quad_start:]
    # w_quad: the calculated weight for the quadratic terms
    w_quad = np.dot(quad_terms, quad_list)

    
    # Total weight
    total_weight = c0 + w_linear + w_quad
    
    return total_weight

def print_structure_constants(struct_constants, event_idx, max_to_print=5):
    # Print the first few structure constants for an event
    wc_names = get_wc_names()
    num_WCs = len(wc_names)
    
    print "\nStructure constants for event {} (showing first {}):".format(event_idx, max_to_print)
    print "Constant term (SM): {}".format(struct_constants[0])
    
    print "\nLinear terms (first few):"
    for i in range(min(max_to_print, num_WCs)):
        print "  {} : {}".format(wc_names[i], struct_constants[1 + i])
    
    if max_to_print < num_WCs:
        print "  ..."

def main():
    parser = argparse.ArgumentParser(description="Test retrieving structure constants and calculating weights")
    parser.add_argument("root_file", help="Path to ROOT file containing structure constants")
    parser.add_argument("--max-events", type=int, default=5, help="Maximum number of events to process")
    parser.add_argument("--verbose", action="store_true", help="Print detailed information")
    parser.add_argument("--wc-name", help="Calculate weight for a specific Wilson coefficient set to 1")
    parser.add_argument("--compare", action="store_true", help="Compare calculated weights with weights from ROOT file")
    args = parser.parse_args()
    
    f = ROOT.TFile.Open(args.root_file)
    if not f or f.IsZombie():
        print "Error: Could not open ROOT file: {}".format(args.root_file)
        return 1
    
    tree = f.Get("AnalysisTree")
    if not tree:
        print "Error: Could not find TTree: AnalysisTree"
        f.Close()
        return 1
    
    # Check if structure_constants branch exists
    branch_names = [b.GetName() for b in tree.GetListOfBranches()]
    if "structure_constants" not in branch_names:
        print "Error: structure_constants branch not found in TTree"
        print "Available branches:", branch_names
        f.Close()
        return 1
    
    # Get number of entries
    nentries = tree.GetEntries()
    print "Found %d events in the ROOT file" % nentries
    
    # Process events
    max_events = min(args.max_events, nentries)
    print "Processing first %d events..." % max_events
    
    # Define Wilson coefficient scenarios
    wc_scenarios = {
        "SM": [0] * 16,
    }
    
    #EDIT SCENARIOS HERE
    # Add scenarios for each Wilson coefficient set to 1
    wc_names = get_wc_names()
    for idx, name in wc_names.items():
        scenario_name = "{} = 1".format(name)
        wc_values = [0] * 16
        wc_values[idx] = 1
        wc_scenarios[scenario_name] = wc_values
    
    # Add more scenarios
    # EDIT SCENARIOS HERE
    wc_scenarios["ctGRe=2"] = [2] + [0] * 15
    wc_scenarios["ctGRe=1,ctGIm=1"] = [1, 1] + [0] * 14
    
    # If a specific WC is requested, only use that one and SM
    if args.wc_name:
        if args.wc_name not in [name for _, name in wc_names.items()]:
            print "Error: Unknown Wilson coefficient name: {}".format(args.wc_name)
            print "Available Wilson coefficients:"
            for _, name in sorted(wc_names.items()):
                print "{}".format(name)
            return 1
        
        # Only use SM and the requested WC
        wc_scenarios = {
            "SM": [0] * 16,
            "{} = 1".format(args.wc_name): create_single_wc_scenario(args.wc_name)
        }
    
    # Process events
    for i in range(max_events):
        tree.GetEntry(i)
        
        # Get structure constants for this event
        struct_constants = np.array(tree.structure_constants)
        
        if args.verbose:
            print_structure_constants(struct_constants, i)
        
        # Calculate weights for different scenarios
        print "\nEvent {} weights:".format(i)
        
        # Get systweights from the ROOT file if we're comparing
        systweights = None
        if args.compare and hasattr(tree, "genInfo"):
            systweights = tree.genInfo.systweights()
        
        for scenario_name, wc_values in sorted(wc_scenarios.items()):
            # Calculate weight using structure constants
            calculated_weight = calculate_weights(struct_constants, wc_values)
            
            # If comparing and we have systweights, get the weight from the ROOT file
            if args.compare and systweights:
                # Get the weight index for this scenario
                weight_idx = get_weight_index(wc_values)
                
                if weight_idx >= 0 and weight_idx < len(systweights):
                    root_weight = systweights[weight_idx]
                    diff = calculated_weight - root_weight
                    diff_percent = 100.0 * diff / root_weight if root_weight != 0 else float('inf')
                    print "  {} : calculated={}, ROOT={} (idx={}), diff={:.2f}%".format(
                        scenario_name, calculated_weight, root_weight, weight_idx, diff_percent)
                else:
                    print "  {} : calculated={}, ROOT=N/A (invalid idx={})".format(
                        scenario_name, calculated_weight, weight_idx)
            else:
                # Just print the calculated weight
                print "  {} : {}".format(scenario_name, calculated_weight)
    
    f.Close()
    return 0

if __name__ == "__main__":
    sys.exit(main()) 