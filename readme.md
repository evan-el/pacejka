# Pacejka

Calculates forces generated at a tire given a .tir file containing parameters of the Pacejka model. Requires slip angles/ratios and vertical force to be provided. The goal is to provide similar functionality to mfeval in Matlab.

## Build
To run the example you will need to provide a .tir file and specify the location in run_example.cpp. Matplot++ is only needed for running the example.

```
mkdir build
cd build
cmake ..
make
./run_example
```

## References
[1] Tyre and Vehicle Dynamics 2nd. ed., Hans B. Pacejka, p. 172-191.

## TODO

Add capability to get tire stiffnesses, not just force.