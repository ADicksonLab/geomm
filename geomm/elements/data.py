# the following method contains portions of the software mdtraj which
# is distributed under the following license
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Peter Eastman, Robert McGibbon
# Contributors: Kyle A. Beauchamp, Matthew Harrigan, Carlos Xavier Hernandez
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
#
# Portions of this code originate from the OpenMM molecular simulation
# toolkit, copyright (c) 2012 Stanford University and Peter Eastman. Those
# portions are distributed under the following terms:
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
##############################################################################


# These van der waals radii are taken from
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
# which references
# A. Bondi (1964). "van der Waals Volumes and Radii". J. Phys. Chem. 68: 441. doi:10.1021/j100785a001 and doi:10.1021/jp8111556.
# M. Mantina et al. (2009). "Consistent van der Waals Radii for the Whole Main Group". J. Phys. Chem. A. 113 (19): 5806--12. doi:10.1021/jp8111556
# Where no van der Waals value is known, a default of 2 angstroms is used.
# However, because certain atoms in biophysical simulations have a high
# chance of being completely ionized, we have decided to give the
# following atoms their ionic radii their ionic radii:
# +2: Be, Mg, Ca, Ba
# +1: Li, Na, K, Cs
# -1: Cl
# These ionic radii are were taken from:
# Shannon, R. D. Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. Acta Crystallographica Section A 32, 751--767 (1976). doi:10.1107/S0567739476001551
# For most atoms, adding electrons usually doesn't change the radius much
# (<10%), while removing them changes it substantially (>50%). Further,
# when atoms like N, S, and P, are positive, they are bound to atoms in such
# a way that would "hide" their radii anyway. We have therefore chosen to just
# use their vdW radii.

ATOMIC_RADII = (('H', 0.12),
                ('He', 0.14),
                ('Li', 0.076),
                ('Be', 0.059),
                ('B', 0.192),
                ('C', 0.17),
                ('N', 0.155),
                ('O', 0.152),
                ('F', 0.147),
                ('Ne', 0.154),
                ('Na', 0.102),
                ('Mg', 0.086),
                ('Al', 0.184),
                ('Si', 0.21),
                ('P', 0.18),
                ('S', 0.18),
                ('Cl', 0.181),
                ('Ar', 0.188),
                ('K', 0.138),
                ('Ca', 0.114),
                ('Sc', 0.211),
                ('Ti', 0.2),
                ('V', 0.2),
                ('Cr', 0.2),
                ('Mn', 0.2),
                ('Fe', 0.2),
                ('Co', 0.2),
                ('Ni', 0.163),
                ('Cu', 0.14),
                ('Zn', 0.139),
                ('Ga', 0.187),
                ('Ge', 0.211),
                ('As', 0.185),
                ('Se', 0.19),
                ('Br', 0.185),
                ('Kr', 0.202),
                ('Rb', 0.303),
                ('Sr', 0.249),
                ('Y', 0.2),
                ('Zr', 0.2),
                ('Nb', 0.2),
                ('Mo', 0.2),
                ('Tc', 0.2),
                ('Ru', 0.2),
                ('Rh', 0.2),
                ('Pd', 0.163),
                ('Ag', 0.172),
                ('Cd', 0.158),
                ('In', 0.193),
                ('Sn', 0.217),
                ('Sb', 0.206),
                ('Te', 0.206),
                ('I', 0.198),
                ('Xe', 0.216),
                ('Cs', 0.167),
                ('Ba', 0.149),
                ('La', 0.2),
                ('Ce', 0.2),
                ('Pr', 0.2),
                ('Nd', 0.2),
                ('Pm', 0.2),
                ('Sm', 0.2),
                ('Eu', 0.2),
                ('Gd', 0.2),
                ('Tb', 0.2),
                ('Dy', 0.2),
                ('Ho', 0.2),
                ('Er', 0.2),
                ('Tm', 0.2),
                ('Yb', 0.2),
                ('Lu', 0.2),
                ('Hf', 0.2),
                ('Ta', 0.2),
                ('W', 0.2),
                ('Re', 0.2),
                ('Os', 0.2),
                ('Ir', 0.2),
                ('Pt', 0.175),
                ('Au', 0.166),
                ('Hg', 0.155),
                ('Tl', 0.196),
                ('Pb', 0.202),
                ('Bi', 0.207),
                ('Po', 0.197),
                ('At', 0.202),
                ('Rn', 0.22),
                ('Fr', 0.348),
                ('Ra', 0.283),
                ('Ac', 0.2),
                ('Th', 0.2),
                ('Pa', 0.2),
                ('U', 0.186),
                ('Np', 0.2),
                ('Pu', 0.2),
                ('Am', 0.2),
                ('Cm', 0.2),
                ('Bk', 0.2),
                ('Cf', 0.2),
                ('Es', 0.2),
                ('Fm', 0.2),
                ('Md', 0.2),
                ('No', 0.2),
                ('Lr', 0.2),
                ('Rf', 0.2),
                ('Db', 0.2),
                ('Sg', 0.2),
                ('Bh', 0.2),
                ('Hs', 0.2),
                ('Mt', 0.2),
                ('Ds', 0.2),
                ('Rg', 0.2),
                ('Cn', 0.2),
                ('Uut', 0.2),
                ('Fl', 0.2),
                ('Uup', 0.2),
                ('Lv', 0.2),
                ('Uus', 0.2),
                ('Uuo', 0.2))


