QuantumLBM
==========

Code for reproducing the results of the 2021 Budinski quantum
advection-diffusion lattice Bolztman paper:
<https://doi.org/10.1007/s11128-021-02996-3>

Files
-----
 * `sriharsha/`
   * `gates.ipynb`: Left shift, Right shift, C1, and C2 gates
   * `AdvectionDiffusion.ipynb`: Circuit for multiple timesteps
   * `Timestep.ipynb`: Circuit for each timestep
   * `classicalLBM.ipynb`: Classical LBM for Advection Diffusion
 * `austin/`
   * `c1.ipynb`: Qiskit code for collision operators `C_1`, `C_2`
   * `c_12.nb`: Mathematica code for collision operators `C_1`, `C_2` generalized to n qubits
   * `cnot-fun-times.nb`: Mathematica code for decomposing
     collision+propagation operators into standard logical quantum gates
   * `L_R.nb`: Mathematica notebook with general decompisition of propagation
     operators `L`, `R`
   * `l_r.ipynb`: Qiskit code showing general decompisition of propagation
     operators `L`, `R`
   * `composition.nb`: Mathematica notebook with final circuit
   * `composition.ipynb`: Qiskit code with final circuit
 * `fatima/`
   * `advection_diffusion_LBM.m`: Classical LBM for advection-diffusion equation using D1Q2. 
      * Set flHeatDiffusion to 1 to solve heat diffusion equation
      * Set flHeatDiffusion to 0 to solve advection-diffusion equation with point source initial condition, compared to Budinski solution
