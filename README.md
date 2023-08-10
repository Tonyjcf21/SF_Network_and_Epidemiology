# Spread of Diseases in SF Networks

This folder contains the code to create scale free networks. The next step would be to modularize the code to make it simpler, 
but for the purpose of my research, it is not necessary, and will be implemented in the future.


- The file **sf_network.cpp** contains the code to generate SF networks.
- The file **sf_network_rewire.cpp** has the same code as the previous mentioned file, but also a section to rewire nodes that are chosen
with a probability by the amount of edges they have.
- The file **sf_network_rewire_SIR2.cpp** contains both previous mentioned functionallities and also a section to spread a disease
based on the SIR epidemiological model, but the spreading needs some adjusting.
- Finally, **sf_network_rewire_SIS2_recoverytime.cpp** and **sf_network_rewire_SIR2_recoverytime.cpp** are both SIS and SIR models that I worked with.

The prints for the output file in the last code are yet to be written but all three codes are working.

