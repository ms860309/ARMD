#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contains constants and unit conversions.
"""
import scipy.constants as spc

###############################################################################

kcal_to_J = 4184
cal_to_J = 4.184
# Speed of light in m/s
C = spc.c
# ħ = h / 2π in J/s
HBAR = spc.hbar
# Conversion factor for momentum in SI to atomic units
P_AU = spc.value("atomic unit of momentum")

# Bohr radius in m
BOHR2M = spc.value("Bohr radius")
# Bohr -> Å conversion factor
BOHR2ANG = BOHR2M * 1e10
# Å -> Bohr conversion factor
ANG2BOHR = 1 / BOHR2ANG
# Hartree to J
AU2J = spc.value("Hartree energy")
# Hartree to kJ / mol
AU2KJPERMOL = AU2J / 1000 * spc.Avogadro
# Hartree to kcal / mol
AU2KCALMOL = AU2KJPERMOL / spc.calorie
# cal to j
CAL2J = AU2KJPERMOL/AU2KCALMOL
# kcal to j
KCAL2J = CAL2J * 1000
# Hartree to eV
AU2EV = spc.value("Hartree energy in eV")
# eV/Å -> Hartree/Bohr
EVANG2AUBOHR = 1 / AU2EV / ANG2BOHR
# fs -> Bohr * sqrt(amu/Hartree)
FS2AU = 0.9682885864793366
# Boltzman constant, (m² kg s⁻² K⁻¹) or just (J / K)
KB = spc.Boltzmann
KBAU = KB / AU2J
# Atomic mass unit to kg
AMU2KG = spc.value("atomic mass constant")
# Atomic unit of time to second
AU2SEC = spc.value("atomic unit of time")
M_E = spc.electron_mass
AMU2AU = AMU2KG / M_E