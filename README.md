QuantumLBM
==========

Code for reproducing the results of the 2021 Budinski quantum
advection-diffusion lattice Bolztman paper:
<https://doi.org/10.1007/s11128-021-02996-3>

Files
-----
 * `sriharsha/`
   * `LRgates.ipynb`: Left and Right shift operators
   * `c1c2.ipynb`: C1 and C2 gates
 * `austin/`
   * `c1.ipynb`: Qiskit code for collision operators `C_1`, `C_2`
   * `cnot-fun-times.nb`: Mathematica code for decomposing
     collision+propagation operators into standard logical quantum gates
   * `L_R.nb`: Mathematica notebook with general decompisition of propagation
     operators `L`, `R`
   * `l_r.ipynb`: Qiskit code showing general decompisition of propagation
     operators `L`, `R`
