# -*- coding: utf-8 -*-

"""
This is a drop-in replacement for "from sympy.physics.units import *"
It loads all the Sympy packages except defintions, which we handle ourselves.
The reason for that is https://github.com/sympy/sympy/issues/14734
"""

from sympy.core.compatibility import string_types
from sympy.physics.units.dimensions import Dimension, DimensionSystem
from sympy.physics.units.unitsystem import UnitSystem
from sympy.physics.units.util import convert_to
from sympy.physics.units.quantities import Quantity

from sympy.physics.units.dimensions import (
    amount_of_substance, acceleration, action,
    capacitance, charge, conductance, current, energy,
    force, frequency, impedance, inductance, length,
    luminous_intensity, magnetic_density,
    magnetic_flux, mass, momentum, power, pressure, temperature, time,
    velocity, voltage, volume,
)

Unit = Quantity

speed = velocity
luminosity = luminous_intensity
magnetic_flux_density = magnetic_density
amount = amount_of_substance

from sympy.physics.units.prefixes import (
    # 10-power based:
    yotta,
    zetta,
    exa,
    peta,
    tera,
    giga,
    mega,
    kilo,
    hecto,
    deca,
    deci,
    centi,
    milli,
    micro,
    nano,
    pico,
    femto,
    atto,
    zepto,
    yocto,
    # 2-power based:
    kibi,
    mebi,
    gibi,
    tebi,
    pebi,
    exbi,
)

from .definitions import (
    percent, percents,
    permille,
    rad, radian, radians,
    deg, degree, degrees,
    sr, steradian, steradians,
    mil, angular_mil, angular_mils,
    m, meter, meters, metre, metres,
    kg, kilogram, kilograms,
    s, second, seconds,
    A, ampere, amperes,
    K, kelvin, kelvins,
    mol, mole, moles,
    cd, candela, candelas,
    g, gram, grams,
    mg, milligram, milligrams,
    ug, microgram, micrograms,
    newton, newtons, N,
    joule, joules, J,
    watt, watts, W,
    pascal, pascals, Pa, pa,
    hertz, hz, Hz,
    coulomb, coulombs, C,
    volt, volts, v, V,
    ohm, ohms,
    siemens, S, mho, mhos,
    farad, farads, F,
    henry, henrys, H,
    tesla, teslas, T,
    weber, webers, Wb, wb,
    optical_power, dioptre, D,
    lux, lx,
    km, kilometer, kilometers,
    dm, decimeter, decimeters,
    cm, centimeter, centimeters,
    mm, millimeter, millimeters,
    um, micrometer, micrometers, micron, microns,
    nm, nanometer, nanometers,
    pm, picometer, picometers,
    ft, foot, feet,
    inch, inches,
    yd, yard, yards,
    mi, mile, miles,
    nmi, nautical_mile, nautical_miles,
    l, liter, liters, litre, litres,
    dl, deciliter, deciliters,
    cl, centiliter, centiliters,
    ml, milliliter, milliliters,
    ms, millisecond, milliseconds,
    us, microsecond, microseconds,
    ns, nanosecond, nanoseconds,
    ps, picosecond, picoseconds,
    minute, minutes,
    h, hour, hours,
    day, days,
    anomalistic_year, anomalistic_years,
    sidereal_year, sidereal_years,
    tropical_year, tropical_years,
    common_year, common_years,
    julian_year, julian_years,
    draconic_year, draconic_years,
    gaussian_year, gaussian_years,
    full_moon_cycle, full_moon_cycles,
    year, years, tropical_year,
    G, gravitational_constant,
    c, speed_of_light,
    Z0,
    hbar,
    planck,
    eV, electronvolt, electronvolts,
    avogadro_number,
    avogadro, avogadro_constant,
    boltzmann, boltzmann_constant,
    R, molar_gas_constant,
    faraday_constant,
    josephson_constant,
    von_klitzing_constant,
    amu, amus, atomic_mass_unit, atomic_mass_constant,
    gee, gees, acceleration_due_to_gravity,
    u0, magnetic_constant,
    e0, electric_constant, vacuum_permittivity,
    Z0, vacuum_impedance,
    coulomb_constant, electric_force_constant,
    atmosphere, atmospheres, atm,
    kPa,
    bar, bars,
    pound, pounds,
    psi,
    dHg0,
    mmHg, torr,
    mmu, mmus, milli_mass_unit,
    quart, quarts,
    ly, lightyear, lightyears,
    au, astronomical_unit, astronomical_units,
    planck_mass,
    planck_time,
    planck_temperature,
    planck_length,
)

