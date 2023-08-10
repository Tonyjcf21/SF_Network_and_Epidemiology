# Spread of Diseases in SF Networks

This folder contains the code to create scale-free networks. The next step would be to modularize the code to make it simpler, 
but for the purpose of my research, it is not necessary and will be implemented in the future.


- The file **sf_network.cpp** contains the code to generate SF networks.
- The file **sf_network_rewire.cpp** has the same code as the previously mentioned file, but also a section to rewire nodes that are chosen
with a probability by the number of edges they have.
- The file **sf_network_rewire_SIS.cpp** contains both previously mentioned functionalities and also a section to spread a disease
based on the SIS epidemiological model.

