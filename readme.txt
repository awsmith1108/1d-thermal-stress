===========================================================
README - 1D Thermal Stress Finite Element Program
===========================================================

This program performs a 1D finite element analysis (FEA) of 
a bar system under combined mechanical loading and thermal 
stresses. It computes nodal displacements, external reactions, 
element strains, and stresses.

-----------------------------------------------------------
1. Input File Format (.dat)
-----------------------------------------------------------

The program reads input from a plain text .dat file that 
defines the nodes, elements, material properties, boundary 
conditions, and temperature distribution.  

The file has three sections:

    1. Header
    2. Node definitions
    3. Element definitions

-----------------------
Header
-----------------------
Format:
    numnod numel
    ref_temp

Where:
    numnod   = number of nodes
    numel    = number of elements
    ref_temp = reference temperature (°C). 
               Thermal strains are assumed zero at this temperature.

-----------------------
Node Definitions
-----------------------
One line per node (numnod total):

    node   xcoord   bcflag   bcvalue

Where:
    node    = node number (system ID)
    xcoord  = node position along x-axis (m)
    bcflag  = boundary condition flag
              0 → prescribed nodal force (N)
              1 → prescribed displacement (m)
    bcvalue = value of prescribed force or displacement

-----------------------
Element Definitions
-----------------------
One line per element (numel total):

    elem   node1   node2   dia   E   temp   alpha

Where:
    elem   = element number (system ID)
    node1  = start node number
    node2  = end node number
    dia    = element diameter (m)
    E      = Young’s modulus (N/m²)
    temp   = element temperature (°C)
    alpha  = coefficient of thermal expansion (1/°C)

-----------------------
Example Input File
-----------------------
    4 3
    20.0
    1  2.0   1   0.0
    2 -1.0   0  -1000.0
    3  0.0   0   2000.0
    4  1.0   0   0.0
    1  1  4  0.01  70e9  10.0  20e-6
    2  2  4  0.02  60e9  50.0  15e-6
    3  3  4  0.03  50e9 100.0  10e-6

This example defines:
    - 4 nodes, 3 elements
    - Reference temperature = 20 °C
    - Node 1: fixed displacement = 0.0
    - Node 2: applied force = -1000 N
    - Node 3: applied force = 2000 N
    - Node 4: applied force = 0 N
    - Elements connected to node 4 with different diameters, 
      stiffnesses, and expansion coefficients

-----------------------------------------------------------
2. Program Overview
-----------------------------------------------------------

The program workflow:

1. Input Handling
   - Reads the .dat file and extracts node/element data.

2. Stiffness Matrix Assembly
   - Constructs global stiffness matrix (K).
   - Computes thermal load vector (G).

3. Boundary Conditions
   - Applies prescribed displacements using a modify-in-place 
     method (modInPlace function).
   - Unknown displacements and reaction forces are solved simultaneously.

4. Solution
   - Solves for nodal displacements (u).
   - Back-substitutes to calculate reaction forces (F).

5. Post-Processing
   - Element strain = displacement difference / element length.
   - Element stress = E * (strain - αΔT).
   - Results are printed: nodal displacements/forces, element 
     strains/stresses.

-----------------------------------------------------------
3. Example Output
-----------------------------------------------------------

    node     displacement (m)   external force (N) 
       1               0.0000        -1.000000e+03 
       2        -6.242119e-04        -1.000000e+03 
       3        -3.615202e-04         2.000000e+03 
       4         3.818914e-04                    0 

    element      strain (m/m)          stress (Pa) 
       1        -3.818914e-04        -1.273240e+07 
       2         5.030516e-04         3.183099e+06 
       3         7.434116e-04        -2.829421e+06 

-----------------------------------------------------------
4. Summary
-----------------------------------------------------------

- Input: geometry, material, boundary, and temperature data (.dat)  
- Output: nodal displacements and forces, element strains and stresses  
- Core Feature: includes thermal stress effects in a 1D finite element bar model  

===========================================================
END OF README
===========================================================
