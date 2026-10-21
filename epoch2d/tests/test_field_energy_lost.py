#!/usr/bin/env python

# Copyright (C) 2009-2019 University of Warwick
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import sdfr
import unittest
from . import SimTest

boundaries = ['x_min', 'x_max', 'y_min', 'y_max']


def readlost(sdffile):
    data = sdfr.read(sdffile, dict=True)
    lost = {}
    for b in boundaries:
        lost[b] = float(data['Field Energy Lost/' + b + ' (J)'].data)
    lost['total'] = float(data['Total Field Energy Lost (J)'].data)
    lost['field'] = float(data['Total Field Energy in Simulation (J)'].data)
    lost['injected'] = \
        float(data['Absorption/Total Laser Energy Injected (J)'].data)
    return lost


class test_field_energy_lost(SimTest):

    def test_initial(self):
        lost = readlost('0000.sdf')
        for b in boundaries + ['total']:
            self.assertEqual(lost[b], 0.0)

    def test_total(self):
        for sdffile in ['0001.sdf', '0002.sdf']:
            lost = readlost(sdffile)
            total = sum(lost[b] for b in boundaries)
            self.assertTrue(np.isclose(lost['total'], total, rtol=1e-12))

    def test_directions(self):
        lost = readlost('0002.sdf')
        # Laser energy enters through x_min and leaves through x_max/y_max
        self.assertLess(lost['x_min'], 0.0)
        self.assertGreater(lost['x_max'], 0.0)
        self.assertGreater(lost['y_max'], 0.0)
        self.assertLess(abs(lost['y_min']), 1e-3 * abs(lost['x_min']))

    def test_injected(self):
        # Compare with the cycle-averaged estimate from the absorption
        # diagnostic, which neglects the obliquity of the pulse
        lost = readlost('0002.sdf')
        self.assertTrue(np.isclose(-lost['x_min'], lost['injected'],
                                   rtol=2e-2))

    def test_conservation(self):
        # After the pulse has left, the energy remaining in the domain plus
        # the energy lost should be zero
        lost = readlost('0002.sdf')
        balance = lost['field'] + lost['total']
        self.assertLess(abs(balance), 1e-4 * abs(lost['x_min']))


if __name__ == '__main__':
    unittest.main()
