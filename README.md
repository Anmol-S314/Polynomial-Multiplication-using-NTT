# Polynomial-Multiplication-using-NTT

The outputs have been tested using the Python implementation of NTT by importing sympy.

This code takes in two arrays of size which is a power of 2 and returns the multiplied product using Iterative Radix2 NTT and Inverse NTT.

```ntt_newfin.c``` is the final C implementation with ```Uno_NTT.ino``` being the Arduino implementation and ```ntt_new1.c``` is used to test it against the Python equivalent.
