# Reading Molcas energy level and g-tensor data. 

I want a python script to read output from molcas, such as: 
ce_alpha_ky3f10.out

The parts of interest are the spin-orbit energy levels: 
--------------------

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
The code allocated initially:           59920 bytes of memory for this run.
MLTP :         =   7
               =   2  2  2  2  2  2  2
The Crystal-Field acting on the ground atomic multiplet of Ln = CE is computed.
CHIT :         =  molar magnetic susceptibility is computed
TINT :      nT = 301
          Tmin =   0.000
          Tmax = 300.000
check_commutation:  The input moment passes all three commutation tests.
AFTER fetch_data_RunFile_all
LOW-LYING SPIN-ORBIT ENERGIES:
ENERGY OF THE SPIN-ORBIT STATE (   1) =         0.00000000000000
ENERGY OF THE SPIN-ORBIT STATE (   2) =         0.00000000000000
ENERGY OF THE SPIN-ORBIT STATE (   3) =       365.53358570950320
ENERGY OF THE SPIN-ORBIT STATE (   4) =       365.53358570950320
ENERGY OF THE SPIN-ORBIT STATE (   5) =       579.24581569762995
ENERGY OF THE SPIN-ORBIT STATE (   6) =       579.24581569762995
ENERGY OF THE SPIN-ORBIT STATE (   7) =      2252.63021709217765
ENERGY OF THE SPIN-ORBIT STATE (   8) =      2252.63021709217765
ENERGY OF THE SPIN-ORBIT STATE (   9) =      2708.39399074633366
ENERGY OF THE SPIN-ORBIT STATE (  10) =      2708.39399074633366
ENERGY OF THE SPIN-ORBIT STATE (  11) =      2917.52501247902865
ENERGY OF THE SPIN-ORBIT STATE (  12) =      2917.52501247902865
ENERGY OF THE SPIN-ORBIT STATE (  13) =      3173.83823161191413
ENERGY OF THE SPIN-ORBIT STATE (  14) =      3173.83823161191413
------------------------
These should be stored an an array: 
energies

And printed out as a column: 
Energy Levels
1 0 
2 0 
... 

I then want the magnetic calcualtions for all level that have been calculated. In this case there are 7 pairs of levels: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET 1 ( effective S =  1/2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The pseudospin is defined in the basis of the following spin-orbit states:
spin-orbit state 1; energy(1) =       0.000 cm-1.
spin-orbit state 2; energy(2) =       0.000 cm-1.
Tunnelling splitting:     0.0000000000 cm-1.

    g TENSOR:
--------------------------------------------------------|
    MAIN VALUES    |             MAIN MAGNETIC AXES     |   x , y , z  -- initial Cartesian axes
-------------------|----|   ----- x ------- y ------- z ---|
Xm, Ym, Zm -- main magnetic axes
 gX =   2.46615127814941 | Xm | -0.69181558287750 -0.72207423391772  0.00000000026082 |
 gY =   2.46615107281067 | Ym |  0.72207423391772 -0.69181558287750 -0.00000000013648 |
 gZ =   0.86115668386291 | Zm |  0.00000000027898  0.00000000009391  1.00000000000000 |
--------------------------------------------------------|
CHECK-SIGN parameter =    0.072770
The sign of the product gX * gY * gZ for multiplet 1: > 0.

--------------------------------------------------
I ned the principal g values gX, gY and gZ, 
and the eigenvectors Xm, Ym, and Zm, as in the printout. 

Again these need to be summarised. 

I then want to use the principle values and eigenvectors to calculate the g-tensor itself, using the code in this file: 
~/calculations/f11/ery2o3_c2/pycf/g_recover.py
and print out the g-tensor. 


I want to add another function to this file to normalize the g tensor by squaring then taking the square root. This converts into a standard form. 

In matlab I did this as follows: 
g_squared = g*g
sqrt_g_squared = g_squared^(1/2)

However, I am not sure of the best square root function to use in python. 

So I need a function that does this calculation and the result shold be printed below the 
 "reconstructed g tensor (G = W * diag(gX,gY,gZ) * W^T):"

 I suggest calling it the "Normalized g tensor"





