---
jupyter:
  jupytext:
    formats: ipynb,src//md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.3.3
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# [QuantumNerd tutorials](https://www.youtube.com/playlist?list=PLGntAYRT8AVmQMyurFoncyOdHljqeGU_R)


Start with QuandumNerd's [intro playlist](https://www.youtube.com/playlist?list=PLGntAYRT8AVki7djmjm4lVV2DImVbno3Z) on Quantum Espresso (QE) before jumping into the projects.

I'm going to try and follow QuantumNerd's projects and bring it into Jupyter [using ASE with QE](https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html)


## [Project 1 H2 molecule](https://www.youtube.com/watch?v=OU9i_zeapCU&list=PLGntAYRT8AVmQMyurFoncyOdHljqeGU_R&index=4)

```python
import os
print(os.getenv('LOCAL_SCRATCH'))
```

```python
from ase import Atoms, Atom
#from ase.calculators.espresso import Espresso
from espresso import Espresso
from ase.units import Bohr
```

```python
pseudopotentials = {'H': 'H_ONCV_PBE-1.0.oncvpsp.upf'}

H2 = Atoms([Atom('H', [0, 0, 0]),
            Atom('H', [1, 0, 0])],
           cell=(20*Bohr, 20*Bohr, 20*Bohr))

H2_calc =  Espresso(label="H2", pseudo_dir = "/Users/matt/QE/SSSP_precision_pseudos/", prefix="H2", 
                    pseudopotentials=pseudopotentials, 
                    calculation="relax", # allows our initial guess of atom positions to change
                    ion_dynamics="bfgs", # relates to the relaxing of the ions
                    ibrav=1,             # this is simple cubic lattice
                    conv_thr=1e-8,       # tells QE when to stop
                    ecutwfc=30.0,        # max energy to use in calculations, bigger this is longer simulation is 
                    kpts=(1, 1, 1)       # number of points in kx,ky,kz space
                   )


H2.calc = H2_calc

print('energy = {0} eV'.format(H2.get_potential_energy()))
```

Final atomic positions cannot be read from the output file by ASE apparently. You have to go and look at the end of `H2.pwo`

```
ATOMIC_POSITIONS (angstrom)
H             0.1189051257        0.0000000000        0.0000000000
H             0.8810948743        0.0000000000        0.0000000000
End final coordinates
```


## [Project: 2 Water molecule](https://www.youtube.com/watch?v=qth17pYTnw4&list=PLGntAYRT8AVmQMyurFoncyOdHljqeGU_R&index=5)

```python
pseudopotentials = {'H': 'H_ONCV_PBE-1.0.oncvpsp.upf',
                    'O': 'O.pbe-n-kjpaw_psl.0.1.UPF'
                   }
# If we don't add some asymmetry in the H's then QE will do detect the symmetry and limit how the molcules can relax
H20 = Atoms([Atom('H', [1, 0.01, 0.03]),
            Atom('H', [-1, 0.01, 0.02]), 
            Atom('O', [0, 0, 0])],
           cell=(10*Bohr, 10*Bohr, 10*Bohr))

H20_calc =  Espresso(label="H20", pseudo_dir = "/Users/matt/QE/SSSP_precision_pseudos/", prefix="H20", 
                    pseudopotentials=pseudopotentials, 
                    calculation="relax",   # allows our initial guess of atom positions to change
                    etot_conv_thr = 1e-5,  # total energy convergence threashold, lower is more accurate (cf "relax")
                    forc_conv_thr = 1e-4,  # force convergence threashold, lower is more accurate (cf "relax")
                    ion_dynamics="bfgs",   # relates to the relaxing of the ions
                    ibrav=1,               # this is simple cubic lattice
                    conv_thr=1e-8,         # tells QE when to stop
                    ecutwfc=25.0,          # max energy to use in calculations, bigger this is longer simulation is 
                    kpts=(1, 1, 1)         # number of points in kx,ky,kz space
                   )

H20.calc = H20_calc

print('energy = {0} eV'.format(H20.get_potential_energy()))
```

```python

```
