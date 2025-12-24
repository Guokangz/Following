#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   constants.py
@Time    :   2025/05/14 10:35:11
@Author  :   George Trenins
@Desc    :   Additional physical/mathematical constants used in all calculations
'''


from __future__ import print_function, division, absolute_import
from dvrlib import units

# unit-system definitions
UNITS = "eVAamu"
u = getattr(units, UNITS)()
# masses
mH = u.str2base('1.008 Da')
muH = mH/2
mH2 = mH*2
mD = mH*2
muD = mD/2
mD2 = mD*2
# σ = 0.96364, uniform scaling factor that accounts for the 
# difference between the stretch frequency of gas-phase H2 #
# calculated with our chosen DFT functional 
freq_scaling = 0.96364
