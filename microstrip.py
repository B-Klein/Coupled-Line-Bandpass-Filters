#!/usr/bin/python2.7 -tt

"""
microstrip.py

This is a port to py of the microstrip.c code written by Claudio Girardi
<claudio.girardi@ieee.org> and Gopal Narayanan 

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
"""

import sys
import math
import cmath

class microstrip:
  er= []			# dielectric constant */
  mur= []			# mag. permeability */
  h= []			# height of substrate */
  ht= []			# height to the top of box */
  t= []			# thickness of top metal */
  sigma= []			# Conductivity of the metal */
  tand= []			# Dielectric Loss Tangent */
  rough= []			# Roughness of top metal */
  f= []			# Frequency of operation */
  w= []			# width of line */
  l= []			# length of line */
  Z0_0= []			# static characteristic impedance */
  Z0= []			# characteristic impedance */
  ang_l= []			# Electrical length in angle */
  er_eff_0= []		# Static effective dielectric constant */
  er_eff= []		# Effective dielectric constant */
  mur_eff= []		# Effective mag. permeability */
  w_eff= []			# Effective width of line */
  atten_dielectric= []	# Loss in dielectric (dB) */
  atten_cond= []		# Loss in conductors (dB) */
  skindepth= []		# Skin depth in mils */
  # private params */
  Z0_h_1= []		# homogeneous stripline impedance */
  def __init__():
    pass

def microstrip_Z0(ms):
  """ /* microstrip_Z0() - compute microstrip static impedance */ """
#   gfloat e_r, h, w, h2, h2h, u, t, t_h;
#   gfloat Z0_h_1, Z0_h_r, Z0;
#   gfloat delta_u_1, delta_u_r, q_inf, q_c, q_t, e_r_eff, e_r_eff_t, q;
  e_r = ms.er;
  h = ms.h;
  w = ms.w;
  t = ms.t;
  h2 = ms.ht;
  h2h = h2 / h;
  u = w / h;
  t_h = t / h;

  # compute normalized width correction for e_r = 1.0 */
  delta_u_1 = delta_u_thickness(u, t_h, 1.0);
  # compute homogeneous stripline impedance */
  Z0_h_1 = Z0_homogeneous(u + delta_u_1);
  # compute normalized width corection */
  delta_u_r = delta_u_thickness(u, t_h, e_r);
  u += delta_u_r;
  # compute homogeneous stripline impedance */
  Z0_h_r = Z0_homogeneous(u);

  # filling factor, with width corrected for thickness */
  q_inf = filling_factor(u, e_r);
  # cover effect */
  q_c = delta_q_cover(h2h);
  # thickness effect */
  q_t = delta_q_thickness(u, t_h);
  # resultant filling factor */
  q = (q_inf - q_t) * q_c;

  # e_r corrected for thickness and non homogeneous material */
  e_r_eff_t = e_r_effective(e_r, q);

  # effective dielectric constant */
  e_r_eff = e_r_eff_t * pow(Z0_h_1 / Z0_h_r, 2.0);

  # characteristic impedance, corrected for thickness, cover */
  #   and non homogeneous material */
  Z0 = Z0_h_r / math.sqrt(e_r_eff_t);

  ms.Z0_h_1 = Z0_h_1;		# will be used in attenuation calculations */
  ms.w_eff = u * h;
  ms.er_eff_0 = e_r_eff;
  ms.Z0_0 = Z0;

def delta_u_thickness(u, t_h, e_r):
  """ /* delta_u_thickness - compute the thickness effect on normalized width */ """
#   gfloat delta_u;
  if t_h > 0.0 :
    # correction for thickness for a homogeneous microstrip */
    delta_u = (t_h / math.pi) * math.log(1.0 + (4.0 * M_E) * pow(math.tanh(math.sqrt(6.517 * u)), 2.0) / t_h);
    # correction for strip on a substrate with relative permettivity e_r */
    delta_u = 0.5 * delta_u * (1.0 + 1.0 / math.cosh(math.sqrt(e_r - 1.0)));
  else:
    delta_u = 0.0;

  return delta_u;

def Z0_homogeneous(u):
  """/* Z0_homogeneous() - compute the impedance for a stripline in a homogeneous medium, without cover effects */ """
#   gfloat f, Z0;
  f = 6.0 + (2.0 * math.pi - 6.0) * math.exp(-pow(30.666 / u, 0.7528));
  Z0 = (377.0 / (2.0 * math.pi)) * math.log(f / u + math.sqrt(1.0 + 4.0 / (u * u)));

  return Z0;

def filling_factor(u, e_r):
  """/* filling_factor() - compute the filling factor for a microstrip without cover and zero conductor thickness */ """
#   gfloat a, b, q_inf;
#   gfloat u2, u3, u4;

  u2 = u * u;
  u3 = u2 * u;
  u4 = u3 * u;

  a = 1.0 + math.log((u4 + u2 / 2704) / (u4 + 0.432)) / 49.0 + math.log(1.0 + u3 / 5929.741) / 18.7;
  b = 0.564 * pow((e_r - 0.9) / (e_r + 3.0), 0.053);

  q_inf = pow(1.0 + 10.0 / u, -a * b);

  return q_inf;

def delta_q_cover(h2h):
  """/* delta_q_cover() - compute the cover effect on filling factor */"""
#   gfloat q_c;

  q_c = math.tanh(1.043 + 0.121 * h2h - 1.164 / h2h);
  return q_c;

def delta_q_thickness(u, t_h):
  """/* delta_q_thickness() - compute the thickness effect on filling factor */"""
  q_t = (2.0 * math.log(2.0) / math.pi) * (t_h / math.sqrt(u));
  return q_t;

def e_r_effective(e_r, q):
  """/* e_r_effective() - compute effective dielectric constant from material e_r and filling factor */ """
  
  e_r_eff = 0.5 * (e_r + 1.0) + 0.5 * q * (e_r - 1.0);
  return e_r_eff;

def microstrip_dispersion(ms):
  """/* microstrip_dispersion() - compute frequency dependent parameters of microstrip */ """

  e_r = ms.er;
  h = ms.h;
  w = ms.w;
  f = ms.f;
  Z0_0 = ms.Z0_0;
  e_r_eff_0 = ms.er_eff_0;
  u = w / h;

  # normalized frequency [GHz * mm] */
  f_n = f * h / 1e06;

  P = e_r_dispersion(u, e_r, f_n);
  # effective dielectric constant corrected for dispersion */
  e_r_eff_f = e_r - (e_r - e_r_eff_0) / (1.0 + P);

  D = Z0_dispersion(u, e_r, e_r_eff_0, e_r_eff_f, f_n);
  Z0_f = Z0_0 * D;

  ms.er_eff = e_r_eff_f;
  ms.Z0 = Z0_f;

def e_r_dispersion(u, e_r, f_n):
  """/* e_r_dispersion() - computes the dispersion correction factor for the effective permeability */"""

  P_1 = 0.27488 + u * (0.6315 + 0.525 / pow(1.0 + 0.0157 * f_n, 20.0)) - 0.065683 * math.exp(-8.7513 * u);
  P_2 = 0.33622 * (1.0 - math.exp(-0.03442 * e_r));
  P_3 = 0.0363 * math.exp(-4.6 * u) * (1.0 - math.exp(-pow(f_n / 38.7, 4.97)));
  P_4 = 1.0 + 2.751 * (1.0 - math.exp(-pow(e_r / 15.916, 8.0)));

  P = P_1 * P_2 * pow((P_3 * P_4 + 0.1844) * f_n, 1.5763);

  return P;

def Z0_dispersion(u, e_r, e_r_eff_0, e_r_eff_f, f_n):
  """/* Z0_dispersion() - computes the dispersion correction factor for the characteristic impedance */"""

  R_1 = 0.03891 * pow(e_r, 1.4);
  R_2 = 0.267 * pow(u, 7.0);
  R_3 = 4.766 * math.exp(-3.228 * pow(u, 0.641));
  R_4 = 0.016 + pow(0.0514 * e_r, 4.524);
  R_5 = pow(f_n / 28.843, 12.0);
  R_6 = 22.2 * pow(u, 1.92);
  R_7 = 1.206 - 0.3144 * math.exp(-R_1) * (1.0 - math.exp(-R_2));
  R_8 = 1.0 + 1.275 * (1.0 - math.exp(-0.004625 * R_3 * pow(e_r, 1.674) * pow(f_n / 18.365, 2.745)));
  tmpf = pow(e_r - 1.0, 6.0);
  R_9 = 5.086 * R_4 * (R_5 / (0.3838 + 0.386 * R_4)) * (math.exp(-R_6) / (1.0 + 1.2992 * R_5)) * (tmpf / (1.0 + 10.0 * tmpf));
  R_10 = 0.00044 * pow(e_r, 2.136) + 0.0184;
  tmpf = pow(f_n / 19.47, 6.0);
  R_11 = tmpf / (1.0 + 0.0962 * tmpf);
  R_12 = 1.0 / (1.0 + 0.00245 * u * u);
  R_13 = 0.9408 * pow(e_r_eff_f, R_8) - 0.9603;
  R_14 = (0.9408 - R_9) * pow(e_r_eff_0, R_8) - 0.9603;
  R_15 = 0.707 * R_10 * pow(f_n / 12.3, 1.097);
  R_16 = 1.0 + 0.0503 * e_r * e_r * R_11 * (1.0 - math.exp(-pow(u / 15.0, 6.0)));
  R_17 = R_7 * (1.0 - 1.1241 * (R_12 / R_16) * math.exp(-0.026 * pow(f_n, 1.15656) - R_15));

  D = pow(R_13 / R_14, R_17);

  return D;


# Define a main() function that prints a little greeting.
def main():
  pass

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
