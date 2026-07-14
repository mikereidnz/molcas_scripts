# open-molcas and helper files installation. 

## Installation of open-molcas

Our patched Molcas is here: https://gitlab.com/chilton-group/open-molcas
### go to where you want to install the repo
cd ~/dev 
git clone https://gitlab.com/chilton-group/open-molcas
git checkout v26.06-patched
### pull in lapack to build ourselves: 
/usr/bin/git submodule update --init /home/users/mfr24/dev/open-molcas/External/lapack
cmake -DCMAKE_INSTALL_PREFIX=[path to install]  -DHDF5=ON ..
mkdir build 
cd build

cmake -DCMAKE_INSTALL_PREFIX=~/open-molcas-26.06  -DHDF5=ON ..
make -j 4
make install

# open-molcas environment 
export PATH=/home/users/mfr24/open-molcas-26.06:$PATH
export MOLCAS=/home/users/mfr24/open-molcas-26.06
export MOLCAS_WORKDIR=/tmp
export MOLCAS_MEM=16000

pymolcas --banner

# on the servers: 
export /scatch/$USER 
 
## Python intallation
### This can become dependency hell... 
cd ~/dev 
python3 -m venv ~/molcas_env
source ~/molcas_env/bin/activate

### The following works, but you end up uninstalling and reinstalling some pacages: 
pip install molcas_suite
pip install --upgrade pip
pip install angmom_suite
pip install env_suite
pip install vasp_suite
pip install phonopy==2.22.0
pip install jax==0.6.2

## in .bashrc: 
source ~/mfr24/dev/molcas_scripts.sh 
This will set up the environment for running pymolcas
to activate the python envionment: 
molcas_activate 

You might consider creating a new file: molcas_scripts_local.sh with your preferences. 

###################################################################################3
## Old versions of instructions. 

## molcas

Convert CIF to POSCAR:  save from VESTA as .vasp

first yttrium
Build QM cluster: env_suite cluster --qm_formula 1 YO6 --cluster_cutoff 30 --poscar Y2O3.vasp
OK two xyz files

could do Y13O6 or Y13O42  
can use 99 and it tells us the  max. 

Make charges file: nano charges.dat (list atomic charges in POSCAR order)32 
32 lines of 3, 48 lines of -2

Set up  charges into “basis set”: env_suite charges ../Y2O3.vasp --from_cluster cluster.xyz --charges charges.dat
writes ENV - charge region only

edit cluster.xyz change "Y " to "Er " atoms (at position 0 0 0) without ENV label are quantum. 
qm_region.xyz is just for visualization and checking. 

Create openMOLCAS file: molcas_suite generate_input cluster.xyz Er1 11 7 6 0 --output er_y203_c2.input --high_S_only --skip_magneto --gateway_extra 'basdir=$CurrDir' --decomp RICD_acCD --kirkwood 1e6 30 4

(molcas_env) l001by:~/calculations/molcas/y2o3/site1$ grep ANO er_y203_c2_basis.xyz 
Er1.ANO-RCC-VTZP       0.0000000       0.0000000       0.0000000
O1924.ANO-RCC-VDZP       1.2735757      -1.2642224       1.4934731
O1927.ANO-RCC-VDZP      -1.2735757      -1.2642224      -1.4934731
O1937.ANO-RCC-VDZP       1.6094387      -0.8132666      -1.3783804
O1939.ANO-RCC-VDZP       1.1584829       1.6187921       1.0425174
O1943.ANO-RCC-VDZP      -1.1584829       1.6187921      -1.0425174
O1957.ANO-RCC-VDZP      -1.6094387      -0.8132666       1.3783804

Change input line:  MLTP= 1; 2  to get more g values.   make that 15;  2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 

can comment out early part of script /*  */ format. 

RUN Molcas:  
pymolcas er_y203_c2.input > er_y203_c2.out 2>&1  &di

LOW-LYING SPIN-ORBIT ENERGIES:
ENERGY OF THE SPIN-ORBIT STATE (   1) =         0.00000000000000
ENERGY OF THE SPIN-ORBIT STATE (   2) =         0.00000000000000
ENERGY OF THE SPIN-ORBIT STATE (   3) =        32.76947477796919
ENERGY OF THE SPIN-ORBIT STATE (   4) =        32.76947477796919
ENERGY OF THE SPIN-ORBIT STATE (   5) =        70.77949043680383
ENERGY OF THE SPIN-ORBIT STATE (   6) =        70.77949043680383
ENERGY OF THE SPIN-ORBIT STATE (   7) =        90.26624155141066
ENERGY OF THE SPIN-ORBIT STATE (   8) =        90.26624155141066
ENERGY OF THE SPIN-ORBIT STATE (   9) =       228.54602400651291
ENERGY OF THE SPIN-ORBIT STATE (  10) =       228.54602400651291
ENERGY OF THE SPIN-ORBIT STATE (  11) =       255.26788213423578
ENERGY OF THE SPIN-ORBIT STATE (  12) =       255.26788213423578
ENERGY OF THE SPIN-ORBIT STATE (  13) =       445.69923412625423
ENERGY OF THE SPIN-ORBIT STATE (  14) =       445.69923412625423
ENERGY OF THE SPIN-ORBIT STATE (  15) =       470.34085351376297
ENERGY OF THE SPIN-ORBIT STATE (  16) =       470.34085351376297


Check active space: molcas_suite orbs 1.rasscf.h5 --index 2

Project CFPs from output: angmom_suite proj --molcas_rassi eq.rassi.h5 --model_space 4I --terms cf=L soc=L,S --quax quax.dat --theta --ion Er3+

angmom_suite proj --molcas_rassi er_y203_c2.rassi.h5 --model_space 4I --terms cf=L soc=L,S --theta --ion Er3+ -o proj.hdf5

use -o proj.out for text file. 

h5ls -dv proj.hdf5/cf/L

----------------

larger cluster did not converge properly... not all 4f... 

-------------------
cd site2 (should have changed the c2 to c3i)
For C3i site try atom 32 to start... 
env_suite cluster --qm_formula 32 YO6 --cluster_cutoff 30 --poscar ../Y2O3.vasp
env_suite charges ../Y2O3.vasp --from_cluster cluster.xyz --charges charges.dat
Edit cluster.xyz as above 
Edit .input as above
run molcas as above
molcas_suite orbs er_y203_c2.rasscf.h5 --index 2
angmom_suite proj --molcas_rassi er_y203_c2.rassi.h5 --model_space 4I --terms cf=L soc=L,S --theta --ion Er3+

That doesn't converge properly, giving one orbital not 4f! 

----------------
Try larger cluster for C3i
cd site-c3i-Y13O6/
env_suite cluster --qm_formula 32 Y13O6 --cluster_cutoff 30 --poscar ../Y2O3.vasp
cp ../site-c3i/charges.dat  .
env_suite charges ../Y2O3.vasp --from_cluster cluster.xyz --charges charges.dat
Edit cluster.xyz to change 
Y           0.0000000       0.0000000       0.0000000
to 
Er          0.0000000       0.0000000       0.0000000
molcas_suite generate_input cluster.xyz Er1 11 7 6 0 --output er_y203_c3i.input --high_S_only --skip_magneto --gateway_extra 'basdir=$CurrDir' --decomp RICD_acCD --kirkwood 1e6 30 4

(molcas_env) l001by:~/calculations/molcas/y2o3/site-c3i-Y13O6$ grep ANO *.xyz
er_y203_c3i_basis.xyz:Y1669.ANO-RCC-VDZ       2.9971724       0.0000000      -2.6519561
er_y203_c3i_basis.xyz:Y1670.ANO-RCC-VDZ       0.0000000       2.6519561      -2.9971724
er_y203_c3i_basis.xyz:Y1673.ANO-RCC-VDZ       2.6519561       2.9971724       0.0000000
er_y203_c3i_basis.xyz:Y1675.ANO-RCC-VDZ      -2.9971724       0.0000000       2.6519561
er_y203_c3i_basis.xyz:Y1676.ANO-RCC-VDZ       0.0000000      -2.6519561       2.9971724
er_y203_c3i_basis.xyz:Y1679.ANO-RCC-VDZ      -2.6519561      -2.9971724       0.0000000
er_y203_c3i_basis.xyz:Y1684.ANO-RCC-VDZ      -2.6519561       2.3067397       0.0000000
er_y203_c3i_basis.xyz:Y1687.ANO-RCC-VDZ       2.3067397       0.0000000       2.6519561
er_y203_c3i_basis.xyz:Y1688.ANO-RCC-VDZ       0.0000000      -2.6519561      -2.3067397
er_y203_c3i_basis.xyz:Y1690.ANO-RCC-VDZ       2.6519561      -2.3067397       0.0000000
er_y203_c3i_basis.xyz:Y1693.ANO-RCC-VDZ      -2.3067397       0.0000000      -2.6519561
er_y203_c3i_basis.xyz:Y1694.ANO-RCC-VDZ       0.0000000       2.6519561       2.3067397
er_y203_c3i_basis.xyz:Er1.ANO-RCC-VTZP       0.0000000       0.0000000       0.0000000
er_y203_c3i_basis.xyz:O1948.ANO-RCC-VDZP      -1.3783804       1.0425174       1.4934731
er_y203_c3i_basis.xyz:O1949.ANO-RCC-VDZP      -1.4934731      -1.3783804      -1.0425174
er_y203_c3i_basis.xyz:O1956.ANO-RCC-VDZP       1.0425174      -1.4934731       1.3783804
er_y203_c3i_basis.xyz:O1960.ANO-RCC-VDZP       1.3783804      -1.0425174      -1.4934731
er_y203_c3i_basis.xyz:O1961.ANO-RCC-VDZP      -1.0425174       1.4934731      -1.3783804
er_y203_c3i_basis.xyz:O1968.ANO-RCC-VDZP       1.4934731       1.3783804       1.0425174

pymolcas er_y203_c3i.input > er_y203_c3i.out 2>&1  &




# Old instructions. 
Download CIF

Convert CIF to POSAR:  vasp_suite convert_cif --filename xxx.cif

Build QM cluster: env_suite cluster --qm_formula 80 YO6 --cluster_cutoff 30 --poscar POSCAR

Make charges file: nano charges.dat (list atomic charges in POSCAR order)

Set up  charges into “basis set”: env_suite charges POSCAR --from_cluster cluster.xyz --charges charges.dat

Create openMOLCAS file: molcas_suite generate_input cluster.xyz Er1 11 7 6 0 --output eq.input --high_S_only --skip_magneto --gateway_extra 'basdir=$CurrDir' --decomp RICD_acCD --Kirkwood 1e6 30 4

RUN Molcas

Check active space: molcas_suite orbs er_y203_c2.rasscf.h5 --index 2

Project CFPs from output: angmom_suite proj --molcas_rassi er_y203_c2.rassi.h5 --model_space 4I --terms cf=L soc=L,S --quax quax.dat --theta --ion Er3+

angmom_suite proj --molcas_rassi er_y203_c2.rassi.h5 --model_space 4I --terms cf=L soc=L,S --theta --ion Er3+
