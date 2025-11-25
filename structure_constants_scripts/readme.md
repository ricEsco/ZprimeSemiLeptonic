1) simple_reweight.py - A script that reads structure constants and calculates new weights for different Wilson coefficient scenarios, with plotting capabilities.
- Calculate weights and create plots
simple_reweight.py your_input_file.root --variable Mass_tt

2) calculate_weights.py - A minimal script that reads structure constants and calculates a new weight for a single Wilson coefficient scenario.
-  Calculate weights for a single scenario
calculate_weights.py your_input_file.root --wc-values "ctGRe=0.5,ctGIm=0.5"


3) calculate_multiple_weights.py - A script that calculates multiple weights at once for different Wilson coefficient scenarios.
- Calculate weights for multiple scenarios defined in the script
calculate_multiple_weights.py your_input_file.root


4) calculate_custom_weights.py - A script that calculates weights for custom Wilson coefficient scenarios.
- Calculate weights and create plots
calculate_custom_weights.py your_input_file.root --wc-values "ctGRe=0.5,ctGIm=0.5"