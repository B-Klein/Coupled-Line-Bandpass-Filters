# Coupled-Line-Bandpass-Filters
Design a coupled line bandpass filter using microstrip lines

usage: bandpassmicrostripfilter.py [-h] [--N N] [--fghz FGHZ] [--r R]
                                   [--bw BW] [--er ER] [--t T] [--mt MT]
                                   [--bh BH] [--sg SG] [--rough ROUGH]
                                   [--lt LT] [--mu MU]

Process environmental variables.

optional arguments:
  -h, --help     show this help message and exit
  --N N          Order of the filter (3)
  --fghz FGHZ    freq in GHz
  --r R          ripple in dBs either 0,0.5 or 3, (0)
  --bw BW        Bandwidth as a propotion delta 0-1, (0.05)
  --er ER        epsilon (4.3)
  --t T          dielectric thickness [mm] (1.6)
  --mt MT        microstrip thickness [mm] (35E-3)
  --bh BH        box height [mm] (1E23)
  --sg SG        sigma (4.1E7)
  --rough ROUGH  roughness (0)
  --lt LT        loss tangent (0)
  --mu MU        relative permability mu (1)
