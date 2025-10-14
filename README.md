# Symbolic–Numerical Modeling and Stability Analysis Platform for DC/AC Power Systems (SyNAPS)

As modern power systems become increasingly complex—integrating diverse inverter-based resources and advanced control schemes—there is a growing need for analysis tools that combine analytical rigor with computational flexibility. Traditional purely numerical simulations often lack transparency. SyNAPS bridges this gap by unifying symbolic computation with numerical simulation, enabling scalable modeling, analysis, and control design for DC and AC power systems. 

It has the following key capabilities:  

- ## Hybrid Symbolic–Numerical Modeling 
Integrates SymPy and SciPy to connect symbolic derivations with numerical evaluation. Symbolic models capture system structure and parameter dependencies, while numerical routines reduce the modeling complexity and enable scalable evaluation and simulation. 

- ## Symbolic Linearization 
Automatically derives Jacobian matrices and linearized state-space models from nonlinear differential algebraic equations (DAEs). This enables transparent analysis of how operating points, grid conditions, and control strategies influence small-signal dynamics. 

- ## Numerical Simulation 
Leverages SciPy for numerical linear algebra, eigenvalue analysis, and time-domain simulation. 

- ## Parametric Stability Characterization 
Defines parameter spaces (e.g., droop gains, PLL gains) and certifies stability using algebraic methods such as Linear Matrix Inequalities (LMIs) and convex optimization. Facilitates control co-design and exploration of stability boundaries. 


# Diagram of symbolic-numerical modeling and stability analysis framework
![Alt_text](diag_software.png)