# -*- coding: utf-8 -*-
#   This file is part of PyMSES.
#
#   PyMSES is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyMSES is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.
"""
:mod:`pymses.utils.constants` --- physical units and constants module
=====================================================================

Use case
--------

    >>> from pymses.utils import constants as C
    >>> l = 100.0 * C.kpc
    >>> t = 320.0 C. Myr
    >>> v = l/t
    >>> print "v = %g km/s"% v.express(C.km/C.s)
    v = 305.56 km/s

"""
from numpy import pi, tan, radians
from .unit import Unit, create_unit
from .physical_type import create_phystype


# unit system = ["m", "kg", "s", "K", "rad", "A", "mol", "cd"]


# --------------------------------------------------------------------------------------------------------------------#
# Dimensionless units
none = create_unit('none',       dims=(0, 0, 0, 0, 0, 0, 0, 0), descr="Unscaled dimensionless unit")
percent = create_unit('percent', base_unit=0.01 * none, descr="One hundredth of unity", latex="%")


# --------------------------------------------------------------------------------------------------------------------#
# Time units
s = create_unit('s',             dims=(0, 0, 1, 0, 0, 0, 0, 0), descr="Second : base time unit")
min = create_unit('min',         base_unit=60. * s,         descr="Minute")
hour = create_unit('hour',       base_unit=3600. * s,       descr="Hour")
day = create_unit('day',         base_unit=24. * hour,      descr="Day")
sid_day = create_unit('sid_day', base_unit=86164.09053 * s, descr="Sidereal day : Earth full rotation time",
                      latex="$T_{sid}$")
year = create_unit('year',       base_unit=365.25 * day,    descr="Year", latex="yr")
kyr = create_unit('kyr',         base_unit=1.0E3 * year,    descr="kyr : millenium")
Myr = create_unit('Myr',         base_unit=1.0E6 * year,    descr="Megayear : million years")
Gyr = create_unit('Gyr',         base_unit=1.0E9 * year,    descr="Gigayear : billion years")


# --------------------------------------------------------------------------------------------------------------------#
# Distance units
m = create_unit('m',               dims=(1, 0, 0, 0, 0, 0, 0, 0), descr="Meter : base length unit")
um = create_unit('um',             base_unit=1.0E-6 * m,         descr="Micron", latex="$\mu m$")
cm = create_unit('cm',             base_unit=1.0E-2 * m,         descr="Centimeter")
nm = create_unit('nm',             base_unit=1.0E-9 * m,         descr="Nanometer")
km = create_unit('km',             base_unit=1000.0 * m,         descr="Kilometer")
Angstrom = create_unit('Angstrom', base_unit=1.0E-10 * m,        descr="Angstrom: 10**-10 m", latex="$\AA$")
au = create_unit('au',             base_unit=1.495978707e8 * km, descr="Astronomical unit")
pc = create_unit('pc',             base_unit=au / tan(radians(1. / 3600.)), descr="Parsec")  # 3.085677e16 * m
kpc = create_unit('kpc',           base_unit=1.0E3 * pc,         descr="Kiloparsec")
Mpc = create_unit('Mpc',           base_unit=1.0E6 * pc,         descr="Megaparsec")
Gpc = create_unit('Gpc',           base_unit=1.0E9 * pc,         descr="Gigaparsec")
Rsun = create_unit('Rsun',         base_unit=6.95508E5 * km,     descr="Solar radius", latex="R$_{\odot}$")


# ------------------------------------------------------------------------------------------------------------------- #
# Freauency units
Hz = create_unit('Hz', base_unit=1.0 / s, descr="Hertz : frequency unit")


# ------------------------------------------------------------------------------------------------------------------- #
# Surface units
barn = create_unit('barn', base_unit=1.0E-28 * m**2, descr="barn: surface unit used in HEP")


# ------------------------------------------------------------------------------------------------------------------- #
# Velocity units
km_s = create_unit('km/s', base_unit=km / s, descr="Kilometers per second", latex="km.s$^{-1}$")

# --------------------------------------------------------------------------------------------------------------------#
# Mass units
kg = create_unit('kg',         coeff=1.0, dims=(0, 1, 0, 0, 0, 0, 0, 0), descr="Kilogram : base mass unit")
g = create_unit('g',           base_unit=1.0E-3 * kg,    descr="Gram")
t = create_unit('t',           base_unit=1.0E3 * kg,     descr="Metric ton")
mH = create_unit('mH',         base_unit=1.66E-24 * g,   descr="Hydrogen atomic mass")
Msun = create_unit('Msun',     base_unit=1.9889E30 * kg, descr="Solar mass", latex="M$_{\odot}$")
Mearth = create_unit('Mearth', base_unit=5.9722E24 * kg, descr="Earth mass", latex="M$_{\oplus}$")


# ------------------------------------------------------------------------------------------------------------------- #
# Angular measurements
rad = create_unit('rad', dims=(0, 0, 0, 0, 1, 0, 0, 0),
                  descr="radian: angular measurement (ratio lengh / radius of an arc")
degree = create_unit('deg',          base_unit=(pi / 180.0) * rad, latex=r"$\textdegree$",
                     descr="degree: angular measurement corresponding to 1/360 of a full rotation")
hourangle = create_unit('hourangle', base_unit=15.0 * degree,
                        descr="hour angle: angular measurement with 24 in a full circle")
arcmin = create_unit('arcmin',       base_unit=degree / 60.0, descr="arc minute: 1/60 of a hour angle", latex="'")
arcsec = create_unit('arcsec',       base_unit=degree/3600.0, descr="arc second: 1/60 of an arcminute", latex="\"")
sr = create_unit('sr',               base_unit=rad ** 2, descr="Steradian: solid angle (SI) unit")

# --------------------------------------------------------------------------------------------------------------------#
# Force units
N = create_unit('N',       base_unit=kg * m / s ** 2, descr="Newton : (SI) force unit")
dyne = create_unit('dyne', base_unit=g * cm / s ** 2, descr="dyne : (CGS) force unit")

# --------------------------------------------------------------------------------------------------------------------#
# Temperature unit
K = create_unit('K', dims=(0, 0, 0, 1, 0, 0, 0, 0), descr="Kelvin : base temperature unit")


# --------------------------------------------------------------------------------------------------------------------#
# Energy units
J = create_unit('J',     base_unit=kg * (m / s)**2,   descr="Joule : (SI) energy unit")
W = create_unit('W',     base_unit=J / s,             descr="Watt")
erg = create_unit('erg', base_unit=g * (cm / s) ** 2, descr="erg : (CGS) energy unit")

# --------------------------------------------------------------------------------------------------------------------#
# Pressure unit
barye = create_unit('barye', base_unit=g / (cm * s**2), descr="Barye: (CGS) pressure unit")
Pa = create_unit('Pa',       base_unit=J * m ** -3,     descr="Pascal: (SI) pressure unit")
hPa = create_unit('hPa',     base_unit=1.0E2 * Pa,      descr="Hectopascal")
kPa = create_unit('kPa',     base_unit=1.0E3 * Pa,      descr="Kilopascal")
bar = create_unit('bar',     base_unit=1.0E5 * Pa,      descr="Bar")
atm = create_unit('atm',     base_unit=1.01325 * bar,   descr="atm: atmospheric pressure (101 3525 Pa)")


# --------------------------------------------------------------------------------------------------------------------#
# Electrical units
A = create_unit('A',     dims=(0, 0, 0, 0, 0, 1, 0, 0), descr="Ampere : electric intensity base unit")
C = create_unit('C',     base_unit=A * s,               descr="Coulomb")
e = create_unit('e',     base_unit=1.602176565e-19 * C, descr="e : elementary electric charge carried by a proton")
V = create_unit('V',     base_unit=J / C,               descr="Volt")
Ohm = create_unit('Ohm', base_unit=V * s / C,           descr="Ohm", latex="$\textohm$")
S = create_unit('S',     base_unit=C / (V * s),         descr="Siemens")
F = create_unit('F',     base_unit=C / V,               descr="Farad")


# --------------------------------------------------------------------------------------------------------------------#
# Magnetical units
T = create_unit('T',           base_unit=V * s / m**2,   descr="Tesla")
Gauss = create_unit('Gauss',   base_unit=1.0E-4 * T,     descr="Gauss")
mGauss = create_unit('mGauss', base_unit=1.0E-3 * Gauss, descr="Milligauss", latex="$mG$")
uGauss = create_unit('uGauss', base_unit=1.0E-6 * Gauss, descr="Microgauss", latex="$\mu G$")
Henry = create_unit('Henry',   base_unit=T * m**2 / A,   descr="Henry",      latex="$H$")


# ------------------------------------------------------------------------------------------------------------------- #
# Amount of substance
mol = create_unit('mol', dims=(0, 0, 0, 0, 0, 0, 1, 0), descr="mole: amount of a chemical substance base unit")


# ------------------------------------------------------------------------------------------------------------------- #
# Illumination
cd = create_unit('cd',     dims=(0, 0, 0, 0, 0, 0, 0, 1), descr="Candela: base luminous intensity unit")
lm = create_unit('lm',     base_unit=cd * sr,             descr="Lumen")
lx = create_unit('lx',     base_unit=lm / m**2,           descr="Lux")
Lsun = create_unit('Lsun', base_unit=3.846E26 * W,        descr="Solar luminosity", latex="L$_{\odot}$")


# ------------------------------------------------------------------------------------------------------------------- #
# Spectral density
Jy = create_unit('Jy', base_unit=1.0E-26 * W / (m**2 * Hz), descr="Jansky")


# --------------------------------------------------------------------------------------------------------------------#
# Composite units
c = create_unit('c',       base_unit=2.99792458E8 * m / s,          descr="Speed of light in vacuum")
ly = create_unit('ly',     base_unit=c * year,                      descr="Light year")
eV = create_unit('eV',     base_unit=e * V,                         descr="electron-Volt")
G = create_unit('G',       base_unit=6.67428e-11 * m ** 3 / kg / (s ** 2), descr="Graviational constant")
kB = create_unit('kB',     base_unit=1.3806504E-23 * J / K,         descr="Boltzmann constant")
H = create_unit('H',       base_unit=70.0 * km / s / Mpc,           descr="Hubble's constant")
rhoc = create_unit('rhoc', base_unit=3.0 * H ** 2 / (8.0 * pi * G), descr="Friedmann's universe critical density",
                   latex="$\rho_{c}$")
H_cc = create_unit('H_cc', base_unit=mH / cm ** 3 / 0.76,           descr="Atoms per cubic centimeter", latex="H/cc")
g_cc = create_unit('g_cc', base_unit=g / cm ** 3,                   descr="Gram per cubic centimeter", latex="g/cc")


# --------------------------------------------------------------------------------------------------------------------#
# Planck constant
h = create_unit('h',       base_unit=6.62606957E-34 * J * s, descr="Planck Constant")
hbar = create_unit('hbar', base_unit=h / (2*pi),             descr="Reduced Planck constant", latex=r"$\bar{h}$")

# --------------------------------------------------------------------------------------------------------------------#
# Stefan-Boltzmann constants
sigmaSB = create_unit('sigmaSB', base_unit=(2 * pi**5 * kB**4) / (15 * h**3 * c**2), descr="Stefan-Boltzmann constant",
                      latex="$\sigma_{SB}$")
a_R = create_unit('a_r',         base_unit=4 * sigmaSB / c, descr="Radiation constant", latex="$a_{R}$")


# ----- Physical type creation ----- #
create_phystype(none, 'dimensionless')
create_phystype(m, 'length')
create_phystype(m**2, 'area')
create_phystype(m**3, 'volume')
create_phystype(s, 'time')
create_phystype(rad, 'angle')
create_phystype(sr, 'solid angle')
create_phystype(m / s, 'velocity')
create_phystype(m / s**2, 'acceleration')
create_phystype(Hz, 'frequency')
create_phystype(g, 'mass')
create_phystype(mol, 'amount of substance')
create_phystype(K, 'temperature')
create_phystype(N, 'force')
create_phystype(J, 'energy')
create_phystype(Pa, 'pressure')
create_phystype(W, 'power')
create_phystype(kg / m**3, 'mass density')
create_phystype(kg / m**2, 'surface density')
create_phystype(m**3 / kg, 'specific volume')
create_phystype(mol / m**3, 'molar volume')
create_phystype(kg * m / s, 'momentum/impulse')
create_phystype(kg * m**2 / s, 'angular momentum')
create_phystype(rad / s, 'angular speed')
create_phystype(rad / s**2, 'angular acceleration')
create_phystype(g / (m * s), 'dynamic viscosity')
create_phystype(m ** 2 / s, 'kinematic viscosity')
create_phystype(m**-1, 'wavenumber')
create_phystype(A, 'electric current')
create_phystype(C, 'electric charge')
create_phystype(V, 'electric potential')
create_phystype(Ohm, 'electric resistance')
create_phystype(S, 'electric conductance')
create_phystype(F, 'electric capacitance')
create_phystype(C * m, 'electric dipole moment')
create_phystype(A / m**2, 'electric current density')
create_phystype(V / m, 'electric field strength')
create_phystype(C / m**2, 'electric flux density')
create_phystype(C / m**3, 'electric charge density')
create_phystype(F / m, 'permittivity')
create_phystype(T * m**2, 'magnetic flux')
create_phystype(T, 'magnetic flux density')
create_phystype(A / m, 'magnetic field strength')
create_phystype(Henry / m, 'electromagnetic field strength')
create_phystype(Henry, 'inductance')
create_phystype(cd, 'luminous intensity')
create_phystype(lm, 'luminous flux')
create_phystype(lx, 'luminous emittence')
create_phystype(W / sr, 'radiant intensity')
create_phystype(cd / m**2, 'luminance')
create_phystype(Jy, 'spectral flux density')

del pi, tan, radians, create_phystype


# ------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------- Rich docstring ---------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #
_mod_doc = """
Built-in constant and unit list
-------------------------------

   * Base units :
      +-----------+-------------------------------------------------------------------+
      |   Name    |                      Description                                  |
      +===========+===================================================================+
"""
_ns = locals()
_sk = list(_ns.keys())
_sk.sort()
_su = []
for _k in _sk:
    _v = _ns[_k]
    if isinstance(_v, Unit) and _v.is_base_unit():
        _su.append(_k)

for _k in _su:
    _v = _ns[_k]
    _mod_doc += "      |%10s | %65s | \n"% (_v._name, _v._descr)
    _mod_doc += "      +-----------+-------------------------------------------------------------------+\n"

_mod_doc += "\n"

_mod_doc += """
   * Derived units :
      +-----------+---------------+---------------------+---------------------------------------------------------""" +\
            """---------------+
      |   Name    |      Value    |         Unit        |                         Description                     """ +\
            """               |
      +===========+===============+=====================+=========================================================""" +\
            """===============+
"""

_su = []
for _k in _sk:
    _v = _ns[_k]
    if isinstance(_v, Unit) and not _v.is_base_unit():
        _su.append(_k)

for _k in _su:
    _v = _ns[_k]
    _mod_doc += "      |%10s |%14g | %19s | %70s | \n"% (_v._name, _v._coeff, _v._decompose_base_units(), _v._descr)
    _mod_doc += "      +-----------+---------------+---------------------+" \
                "------------------------------------------------------------------------+\n"

_mod_doc += "\n"

__doc__ += _mod_doc

del _mod_doc, _k, _v, _su, _sk, _ns
# ------------------------------------------------------------------------------------------------------------------- #