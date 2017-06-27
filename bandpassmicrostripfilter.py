#!/usr/bin/python2.7 -tt

"""
bandpassmicrostripfilter.py

Copyright (C) 2017 Benjamin Klein <Benjamin.a.Klein@protonmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA, 02111-1307, USA

Creates a bandpass filter based on Pozar's Chapter 8 filter design method, and the
microstrip impedence calculations done in the open source program Transcalc.

"""

""" TODO:
  Don't hard code the filter coefficients ripple should be a variable.
"""

import sys
import math
import cmath
import c_microstrip as cms
import argparse

# Zo = 50
c = 299792458
theta = math.pi/2

""" c and L parameters 
from Pozar : TABLE 8.4-5 Element Values for Equal-Ripple Low-Pass Filter Prototypes (g0 = 1, wc =
1, N = 1 to 10, 0, 0.5 dB and 3.0 dB ripple) """
#0 dB ripple
db0ripple = {1:[2.0000, 1.0000]}
db0ripple[2] = [1.5774, 0.4226, 1.0000]
db0ripple[3] = [1.2550, 0.5528, 0.1922, 1.0000]
db0ripple[4] = [1.0598, 0.5116, 0.3181, 0.1104, 1.0000]
db0ripple[5] = [0.9303, 0.4577, 0.3312, 0.2090, 0.0718, 1.0000]
db0ripple[6] = [0.8377, 0.4116, 0.3158, 0.2364, 0.1480, 0.0505, 1.0000]
db0ripple[7] = [0.7677, 0.3744, 0.2944, 0.2378, 0.1778, 0.1104, 0.0375, 1.0000]
db0ripple[8] = [0.7125, 0.3446, 0.2735, 0.2297, 0.1867, 0.1387, 0.0855, 0.0289, 1.0000]
db0ripple[9] = [0.6678, 0.3203, 0.2547, 0.2184, 0.1859, 0.1506, 0.1111, 0.0682, 0.0230, 1.0000]
db0ripple[10] = [0.6305, 0.3002, 0.2384, 0.2066, 0.1808, 0.1539, 0.1240, 0.0911, 0.0557, 0.0187, 1.0000]
#0.5 dB ripple
db0p5ripple = {1:[0.6986, 1.0000]}
db0p5ripple[2] = [1.4029,0.7071,1.9841]
db0p5ripple[3] = [1.5963,1.0967,1.5963,1.0000]
db0p5ripple[4] = [1.6703,1.1926,2.3661,0.8419,1.9841]
db0p5ripple[5] = [1.7058,1.2296,2.5408,1.2296,1.7058,1.0000]
db0p5ripple[6] = [1.7254,1.2479,2.6064,1.3137,2.4758,0.8696,1.9841]
db0p5ripple[7] = [1.7372,1.2583,2.6381,1.3444,2.6381,1.2583,1.7372,1.0000]
db0p5ripple[8] = [1.7451,1.2647,2.6564,1.3590,2.6964,1.3389,2.5093,0.8796,1.9841]
db0p5ripple[9] = [1.7504,1.2690,2.6678,1.3673,2.7239,1.3673,2.6678,1.2690,1.7504,1.0000]
db0p5ripple[10] = [1.7543,1.2721,2.6754,1.3725,2.7392,1.3806,2.7231,1.3485,2.5239,0.8842,1.9841]
#3 dB ripple
db3ripple = {1:[1.9953, 1.0000]}
db3ripple[2]=[3.1013, 0.5339, 5.8095]
db3ripple[3]=[3.3487, 0.7117, 3.3487, 1.0000]
db3ripple[4]=[3.4389, 0.7483, 4.3471, 0.5920, 5.8095]
db3ripple[5]=[3.4817, 0.7618, 4.5381, 0.7618, 3.4817, 1.0000]
db3ripple[6]=[3.5045, 0.7685, 4.6061, 0.7929, 4.4641, 0.6033, 5.8095]
db3ripple[7]=[3.5182, 0.7723, 4.6386, 0.8039, 4.6386, 0.7723, 3.5182, 1.0000]
db3ripple[8]=[3.5277, 0.7745, 4.6575, 0.8089, 4.6990, 0.8018, 4.4990, 0.6073, 5.8095]
db3ripple[9]=[3.5340, 0.7760, 4.6692, 0.8118, 4.7272, 0.8118, 4.6692, 0.7760, 3.5340, 1.0000]
db3ripple[10]=[3.5384, 0.7771, 4.6768, 0.8136, 4.7425, 0.8164, 4.7260, 0.8051, 4.5142, 0.6091, 5.8095]

def getgconstants(n,N,ripple):

  global db0ripple
  global db0p5ripple
  global db3ripple
    
  if N > 10: 
    sys.stderr.write('N must be 10 or less') 
    sys.exit(1)

  if ripple == 0:
    g = db0ripple[N][n-1]
  elif ripple == 0.5:
    g = db0p5ripple[N][n-1]
  elif ripple == 3:
    g = db3ripple[N][n-1]
  else:
    sys.stderr.write('only ripples allowed: 0,0.5,3')
    sys.exit(1)

  return g

def getZJconstants(n,N,ripple,deltabw):
  """ implements Pozar 8.121a-c """
  
  if n == 1:
    ZoJn = math.sqrt((math.pi*deltabw)/(2*getgconstants(n,N,ripple)))
    return ZoJn
  elif n == (N+1):
    ZoJn = math.sqrt((math.pi*deltabw)/(2*getgconstants(N,N,ripple)*getgconstants(N+1,N,ripple)))
    return ZoJn
  else:
    ZoJn = (math.pi*deltabw)/(2*math.sqrt(getgconstants(n-1,N,ripple)*getgconstants(n,N,ripple)))
    return ZoJn
  return None

def createfilter(deltabw,N,ripple):

  global Zo
  global theta
  ZJ = []; Zoe = []; Zoo = []; cosB=[]; B=[]; L=[]
  for n in range(1,N+2):
    zj = getZJconstants(n,N,ripple,deltabw)
    zoe = Zo*(1 + zj + pow((zj),2))
    zoo = Zo*(1 - zj + pow((zj),2))
    ZJ.append(zj)
    Zoe.append(zoe)
    Zoo.append(zoo)
#     cosb = (zj + 1/zj)*math.sin(theta)*math.cos(theta)
#     B.append(math.acos(cosb))
#     cosB.append(cosb)
#     L.append(theta/(math.acos(cosb)))
# 
#   print L
  return Zoo,Zoe


# Define a main() function that prints a little greeting.
def main():
  """ Runs the code to synthasise a bandpass stepped filter"""

  global c
  global theta
  global Zo
  
  parser = argparse.ArgumentParser(description='Process environmental variables.')
  parser.add_argument('--N', default=3, type=int, help='Order of the filter (3)')
  parser.add_argument('--fghz', default=1, type=float, help='freq in GHz')
  parser.add_argument('--r', default=0, type=float, help='ripple in dBs either 0,0.5 or 3, (0)')
  parser.add_argument('--bw', default=0.05, type=float, help='Bandwidth as a propotion delta 0-1, (0.05)')
  parser.add_argument('--er', default=4.3, type=float, help='epsilon (4.3)')
  parser.add_argument('--t', default=1.6, type=float, help='dielectric thickness [mm] (1.6)')
  parser.add_argument('--mt', default=35E-3, type=float, help='microstrip thickness [mm] (35E-3)')
  parser.add_argument('--bh', default=1E23, type=float, help='box height [mm] (1E23)')
  parser.add_argument('--sg', default=4.1E7, type=float, help='sigma (4.1E7)')
  parser.add_argument('--rough', default=0, type=float, help='roughness (0)')
  parser.add_argument('--lt', default=0, type=float, help='loss tangent (0)')
  parser.add_argument('--mu', default=1, type=float, help='relative permability mu (1)')
  parser.add_argument('--z', default=50, type=float, help='characteristic impedance (50)')
  args = parser.parse_args()
  # unpack the parser
  deltabw = args.bw
  N = args.N
  ripple = args.r
  freqGHz = args.fghz
  Zo = args.z

  fr4 = cms.cmicrostrip
  fr4.er = args.er
  fr4.h = args.t*1E-3 # 1.6E-3 # mm2mil(1.6)
  fr4.t = args.mt*1E-3 # 35E-6 #umm2mil(35)
  fr4.mur = args.mu # 1
#   fr4.ht = cms.mil2mm(1E+20)*1E-3
  fr4.ht = args.bh*1E-3 #1E+20
  fr4.f = freqGHz*1E9
  fr4.sigma = args.sg #4.1e7 # 4.1E2
  fr4.rough = args.rough # 0
  fr4.tand = args.lt # 0
  fr4.l = cms.mil2mm(1000)*1E3
  wmm = smm = []
  
  # Let's calculate the dimensions of the microstrip lines
  lmdmm = (1E3*c)/(freqGHz*1E9)
  [Zoo,Zoe]=createfilter(deltabw,N,ripple)
# 
  print 'lmd/4 [mm]:',lmdmm/4
  B = 1/math.sqrt(fr4.er)
  print 'theta=B*L [mm]:',B*lmdmm/4
  print 'using 45 deg squaresize [mm]:',1/math.sqrt(2)*B*lmdmm/4

  for n in range(0,len(Zoe)):
    print '\nn =',n
    print 'Zoe =',Zoe[n],'Zoo =',Zoo[n]
    fr4.Z0e = Zoe[n]
    fr4.Z0o = Zoo[n]
    try:
      cms.synthesize_c_microstrip(fr4)
    except:
      print 'filter is not physically realisable, maybe broaden the bandwidth or reduce the order?'
      sys.exit(1)
    wmm.append(fr4.w*1E3)
    smm.append(fr4.s*1E3)
    print 'w [mm] =',fr4.w*1E3, 's [mm] =',fr4.s*1E3

#   print 'er_eff_e',fr4.er_eff_e
#   print 'er_eff_o',fr4.er_eff_o
#   print 'Conductor losses even [dB]',fr4.atten_cond_e*1E-6
#   print 'Conductor losses odd [dB]',fr4.atten_cond_o*1E-6
#   print 'Dielectric losses even [dB]',fr4.atten_dielectric_e*1E-6
#   print 'Dielectric losses odd [dB]',fr4.atten_dielectric_o*1E-6
#   print 'Skin depth um', fr4.skindepth*1E6


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
