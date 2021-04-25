# -*- coding: utf-8 -*-
"""

"""

from ase.calculators.vasp import Vasp
from ase.io import write
from ase.io import read
from ase.constraints import constrained_indices
from ase.constraints import FixAtoms
import numpy as np
import os
import sys
import time
from copy import deepcopy
import subprocess
import random

proc = subprocess.Popen('whoami'.split(), stdout=subprocess.PIPE)
username = (proc.communicate())[0].decode().strip()

def submit_job(job_path, jobname):
    """
    Submit a job to cluster with PBS scheduler
    """
    # Submit a job
    proc = subprocess.Popen('qsub -N {} {}'.format(jobname, job_path),
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    # Get job ID
    stdout = [line.decode().rstrip() for line in proc.stdout]
    job_id = int(stdout[0].split('.')[0])
    proc.wait()

    # Wait till job finish
    finish = False
    while not finish:
        proc = subprocess.Popen('qstat -u {un}'.format(un=username),
                                shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout = [line.decode().rstrip() for line in proc.stdout]
        running_jobs_ids = list()
        for i, status_line in enumerate(stdout):
            if i > 4:
                job_id_ = int(status_line.strip().split('.')[0])
                running_jobs_ids.append(job_id_)
        if job_id not in running_jobs_ids:
            finish = True
        time.sleep(30)
    # Clear job files
    subprocess.Popen(['rm {jobname}.o{ID} {jobname}.e{ID}'.format(jobname=jobname, ID=job_id)], shell=True)


def run(jobname):
    """
    Run a scf job and a Chargemol program
    """
    # Get path of working directory
    path = os.getcwd()
    # Get job file path
    job_path = os.path.join(path, 'job.sh')
    # Submit a job to do a scf calculation
    submit_job(job_path, 'vasp_'+jobname)



# CO coordinate 
f = open('./POSCAR','r+')
line_to_replace = 17
data = f.readlines()
Ce_coord = data[14].split()
C_z = random.uniform(float(Ce_coord[2])+1.5, float(Ce_coord[2])+2.5)
C_y = random.uniform(float(Ce_coord[1]), float(Ce_coord[1])+2)
C_x = random.uniform(float(Ce_coord[0]), 3.8775703907)

dist_C_Ce = ((float(Ce_coord[2])-C_z)**2+(float(Ce_coord[1])-C_y)**2+(float(Ce_coord[0])-C_x)**2)**0.5
while dist_C_Ce < 2.5 :
    C_z = random.uniform(float(Ce_coord[2])+1.5, float(Ce_coord[2])+2.5)
    C_y = random.uniform(float(Ce_coord[1]), float(Ce_coord[1])+2)
    C_x = random.uniform(float(Ce_coord[0]), 3.8775703907)
    dist = ((float(Ce_coord[2])-C_z)**2+(float(Ce_coord[1])-C_y)**2+(float(Ce_coord[0])-C_x)**2)**0.5

        
O_z = random.uniform(C_z, C_z+1.5)
O_y = random.uniform(C_y, C_y+1.5)
O_x = random.uniform(C_x, C_x+1.5)

dist = ((O_z-C_z)**2+(O_y-C_y)**2+(O_x-C_x)**2)**0.5
while dist < 1.12 or dist > 1.15:
    O_z = random.uniform(C_z, C_z+1.5)
    O_y = random.uniform(C_y, C_y+1.5)
    O_x = random.uniform(C_x, C_x+1.5)
    dist = ((O_z-C_z)**2+(O_y-C_y)**2+(O_x-C_x)**2)**0.5

f.write('\n')
f.write(str(O_x)+' '+str(O_y)+' '+str(O_z)+' '+'T'+' '+'T'+' '+'T')
f.write('\n')
f.write(str(C_x)+' '+str(C_y)+' '+str(C_z)+' '+'T'+' '+'T'+' '+'T')
f.close()

# Read file
slab = read('./POSCAR')
slab_save = deepcopy(slab)
isconstrained = constrained_indices(slab)
All_Atoms = np.array(range(0,len(slab)))
Free_Atoms = [atom for atom in All_Atoms if atom not in isconstrained]

# Run the scf job at stationary potition
calc = Vasp(gga='PE', xc='PBE', gamma=True, kpts=(6,6,1), lwave=False, laechg=False, encut=450, ismear=0, sigma=0.05, nsw=300, prec='Accurate', nelmin=4, nelm=400, ediff=1e-8, ediffg=-0.01, ispin=1, isif=2, ibrion=2, npar=4, lreal='Auto', idipol=3, ldau=True, ldautype=2, ldaul=[3,-1,-1], ldauu=[5,0,0], ldauj=[0,0,0], ldauprint=2, lmaxmix=6, ivdw=11)
calc.initialize(slab)
calc.write_kpoints()
#calc.write_potcar()
calc.write_incar(slab)
write('POSCAR', calc.atoms_sorted)
run(jobname='110_surf')











