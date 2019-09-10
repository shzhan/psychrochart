# import os
# os.chdir("./psychrochart")
# sys.path.append(os.getcwd())

# -*- coding: utf-8 -*-
"""A library to make psychrometric charts and overlay information in them."""
from typing import Any
from math import log, exp, atan, sqrt
from scipy.optimize import fsolve

from psychrochart.util import iter_solver


PRESSURE_STD_ATM_PSIA = 14.696  #CHECKED, eq 3
DELTA_TEMP_F_TO_RANKINE = 459.67  #CHANGED 

# Eq. (1) 2009 ASHRAE Handbook—Fundamentals (IP)
GAS_CONSTANT_R_DA = 53.350 # ft/lb_f/lb_da R
# Eq. (2) 2009 ASHRAE Handbook—Fundamentals (SI)
GAS_CONSTANT_R_W = 85.780  # = 8314.472/18.015268 = J/(kgw·K)


def pressure_by_altitude(altitude_ft: float) -> float:#CHANGED

    """Obtain the standard pressure for a certain sea level.

    Pressure by altitude, eq. (3) 2017 ASHRAE Handbook—Fundamentals (IP).
    """
    return PRESSURE_STD_ATM_PSIA * (1 - 6.8754**(-6) * altitude_ft) ** 5.2559  # in psia

def water_vapor_pressure(
        w_lb_lba: float, p_atm_psia: float=PRESSURE_STD_ATM_PSIA) -> float:  # Changed
    """Obtain the water vapor pressure from the humidity ratio (w_kg_kga).


    from eq 20 of 2017 ASHRAE Handbook—Fundamentals (IP).
    """
    humid_ratio_vap_pres = .621945
    p_vapor_psia = p_atm_psia * w_lb_lba / (humid_ratio_vap_pres + w_lb_lba)
    return p_vapor_psia


def humidity_ratio(
        p_vapor_psia: float, p_atm_psia: float=PRESSURE_STD_ATM_PSIA) -> float:
    """Obtain the humidity ratio from the water vapor pressure.

     eq (20)
    2017 ASHRAE Handbook—Fundamentals (IP).
    """
    humid_ratio_vap_pres = .621945
    w_lb_lba = humid_ratio_vap_pres * p_vapor_psia / (p_atm_psia - p_vapor_psia)
    return w_lb_lba


# TODO prec revision:
def humidity_ratio_from_temps(  #CHANGED
        dry_bulb_temp_F: float, wet_bulb_temp_F: float,
        p_atm_psia: float=PRESSURE_STD_ATM_PSIA) -> float:
    """Obtain the specific humidity from the dry and wet bulb temperatures.

    kg water vapor / kg dry air, eqs (33) and (35)
    2017 ASHRAE Handbook—Fundamentals (IP).
    """
    w_sat_wet_bulb = humidity_ratio(
        saturation_pressure_water_vapor(wet_bulb_temp_F), p_atm_psia)
    delta_t = dry_bulb_temp_F - wet_bulb_temp_F
    # assert (delta_t >= 0)
    factor_delta = 0.24 * delta_t
    if dry_bulb_temp_F > 32:
        num_1 = (1093 - 0.556 * wet_bulb_temp_F) * w_sat_wet_bulb
        denom_1 = 1093 + 0.444 * dry_bulb_temp_F -  wet_bulb_temp_F
        # print(num_1 - factor_delta)
        w_lb_lba = (num_1 - factor_delta) / denom_1
    else:
        num_2 = (1220 - 0.04 * wet_bulb_temp_F) * w_sat_wet_bulb
        denom_2 = 1220 + 0.444 * dry_bulb_temp_F - 0.48 * wet_bulb_temp_F
        w_lb_lba = (num_2 - factor_delta) / denom_2
    return w_lb_lba


def relative_humidity_from_temps(   #CHANGED  #CHECKED
        dry_bulb_temp_F: float, wet_bulb_temp_F: float,
        p_atm_psia: float=PRESSURE_STD_ATM_PSIA) -> float:
    """Obtain the relative humidity from the dry and wet bulb temperatures.

    Ratio of the mole fraction of water vapor x_w in a given moist
    air sample to the mole fraction xws in an air sample
    saturated at the same temperature and pressure.

    Eq (24) 2009 ASHRAE Handbook—Fundamentals (IP).
    """
    w_lb_lba = humidity_ratio_from_temps(
        dry_bulb_temp_F, wet_bulb_temp_F, p_atm_psia)
    p_w_vapor = water_vapor_pressure(w_lb_lba, p_atm_psia)
    p_sat = saturation_pressure_water_vapor(dry_bulb_temp_F)
    return p_w_vapor / p_sat
    # return min(1., p_w_vapor / p_sat)


def specific_volume(  #CHANGED
        dry_temp_F: float, p_vapor_psia: float,
        p_atm_psia: float=PRESSURE_STD_ATM_PSIA) -> float:
    """Obtain the specific volume v of a moist air mixture.

    m3 / kg dry air, eq. (28) 2009 ASHRAE Handbook—Fundamentals (SI).
    """
    w_lb_lba = humidity_ratio(p_vapor_psia, p_atm_psia)
    specific_vol_ft3_lba = (GAS_CONSTANT_R_DA
                           * (dry_temp_F + DELTA_TEMP_F_TO_RANKINE)
                           * (1 + 1.607858 * w_lb_lba) / p_atm_psia)
    return specific_vol_ft3_lba


def dry_temperature_for_specific_volume_of_moist_air(  ##CHANGED  # check why is it in inHg
        w_lb_lba: float, specific_vol: float,
        p_atm_psia: float=PRESSURE_STD_ATM_PSIA) -> float:
    """Solve the dry bulb temp from humidity ratio and specific volume.

    ºC. Derived from eq. (26), 2017 ASHRAE Handbook—Fundamentals (IP).
    """
    # p_atm_inHg = 2.0360206576012 * p_atm_psia
    # return (1349578.6615419 * p_atm_inHg * specific_vol
    #         / (803929. * w_lb_lba + 500000.) - DELTA_TEMP_F_TO_RANKINE)
    return  specific_vol * p_atm_psia / (GAS_CONSTANT_R_DA * (1 + 1.607858 * w_lb_lba)) - DELTA_TEMP_F_TO_RANKINE


def density_water_vapor(p_vapor_psia: float, dry_bulb_temp_F: float) -> float:
    """Density of water vapor."""
    return (p_vapor_psia
            / ((dry_bulb_temp_F + DELTA_TEMP_F_TO_RANKINE) * GAS_CONSTANT_R_W))


def saturation_pressure_water_vapor(dry_temp_F: float,
                                    mode: int=1,
                                    logger: Any=None) -> float:   #CHANGED
    """Saturation pressure of water vapor (kPa) from dry temperature.

    3 approximations:
      - mode 1: ASHRAE formulation
      - mode 2: Simpler, values for T > 0 / T < 0, but same speed as 1
      - mode 3: More simpler, near 2x vs mode 1.
    """
    def _handle_error(exception):  # pragma: no cover
        msg = 'OverflowError: {} with T={:.2f} °C, mode={}. ' \
              'Changing to ASHRAE eqs'.format(exception, dry_temp_F, mode)
        if logger is not None:
            logger.error(msg)  # pragma: no cover
        else:
            print(msg)
        return saturation_pressure_water_vapor(dry_temp_F, mode=1)

    if mode == 1:  # 2009 ASHRAE Handbook—Fundamentals (IP)
        abs_temp = dry_temp_F + DELTA_TEMP_F_TO_RANKINE
        if dry_temp_F > 32:  # Eq (6) 2009 ASHRAE Handbook—Fundamentals (IP)
            c1 = -1.0440397E4
            c2 = -1.1294650E1
            c3 = -2.7022355E-2
            c4 = 1.2890360E-5
            c5 = -2.4780681E-9
            c6 = 6.5459673
            ln_p_ws_psia = (c1 / abs_temp + c2 + c3 * abs_temp
                          + c4 * abs_temp ** 2 + c5 * abs_temp ** 3
                          + c6 * log(abs_temp))
        else:  # Eq (5) 2009 ASHRAE Handbook—Fundamentals (SI)
            c7 = -1.0214165E4
            c8 =  -4.8932428E+00
            c9 = -5.3765794E-03
            c10 = 1.9202377E-07
            c11 = 3.5575832E-10
            c12 = -9.0344688E-14
            c13 = 4.1635019E+00
            # print(abs_temp)
            ln_p_ws_psia = (c7 / abs_temp + c8 + c9 * abs_temp
                          + c10 * abs_temp ** 2 + c11 * abs_temp ** 3
                          + c12 * abs_temp ** 4 + c13 * log(abs_temp))
        return exp(ln_p_ws_psia)
    elif mode == 2:  #Tetens's formula
        dry_temp_c = (dry_temp_F - 32) * 5 / 9
        if dry_temp_F > 32:
            return 0.1450377377*(610.5 * exp(
                17.269 * dry_temp_c / (237.3 + dry_temp_c)) / 1000.)
        else:
            try:
                return 0.1450377377*(610.5 * exp(
                    21.875 * dry_temp_c / (265.5 + dry_temp_c)) / 1000.)
            except OverflowError as exc:  # pragma: no cover
                return _handle_error(exc)
    else:  # mode = 3
        dry_temp_c = (dry_temp_F - 32) * 5 / 9
        try:
            return 0.1450377377 * (19314560 * 10. ** (-1779.75 / (237.3 + dry_temp_c)))
        except OverflowError as exc:  # pragma: no cover
            return _handle_error(exc)


def enthalpy_moist_air(  #CHANGED
        dry_temp_F: float, p_vapor_psia: float,
        p_atm_psia: float=PRESSURE_STD_ATM_PSIA) -> float:
    """Moist air specific enthalpy.

    KJ / kg. Eqs. (32), (30), (31) 2009 ASHRAE Handbook—Fundamentals (SI).
    """
    w_lb_lba = humidity_ratio(p_vapor_psia, p_atm_psia)

    c_pa = 0.240  # kJ/(kg·ºC), sensible
    c_pv = 0.444  # kJ/(kg·ºC), latent
    # Evaporation heat (using constant value. For water at 0ºC, it's 2501)
    # constant_water_evaporation_heat = 2503  # KJ / kg
    constant_water_evaporation_heat = 1061  # KJ / kg

    h_a = c_pa * dry_temp_F
    h_v = w_lb_lba * (constant_water_evaporation_heat + c_pv * dry_temp_F)
    h_m = h_a + h_v
    return h_m


def dry_temperature_for_enthalpy_of_moist_air(   #CHANGED   #CHECKED
        w_lb_lba: float, enthalpy: float) -> float:
    """Solve the dry bulb temp from humidity ratio and specific enthalpy.

    ºC. Derived from eq. (30), 2017 ASHRAE Handbook—Fundamentals (SI).
    """
    return 500. * (enthalpy - 1061 * w_lb_lba) / (222.0 * w_lb_lba + 120.0)


def dew_point_temperature(p_w_psia: float) -> float:
    """Dew point temperature.

    Eqs. (37) and (38) Peppers 1988
    2017 ASHRAE Handbook—Fundamentals (ip).
    """
    try:
        alpha = log(p_w_psia)
    except ValueError:  # pragma: no cover
        raise AssertionError("Bad water vapor pressure! ({} kPa)"
                             .format(p_w_psia))
    c14 = 100.45
    c15 = 33.193
    c16 = 2.319
    c17 =  0.17074
    c18 = 1.2063
    dew_point = (c14 + c15 * alpha + c16 * alpha ** 2 + c17 * alpha ** 3
                 + c18 * p_w_psia ** .1984)
    if dew_point < 32:
        # print('BAD dew_point: ', dew_point)
        dew_point = 90.12 + 26.142 * alpha + .8927 * alpha ** 2
        # print('GOOD dew_point: ', dew_point)
    return dew_point


def wet_bulb_temperature(  #CHANGED
        dry_temp_F: float, relative_humid: float,
        p_atm_psia: float=PRESSURE_STD_ATM_PSIA,
        num_iters_max=10000000, precision=0.1) -> float:
    """Wet bulb temperature.

    Requires trial-and-error or numerical solution method
    applying eqs. (23) and (35) or (37)
    2009 ASHRAE Handbook—Fundamentals (SI).
    """

    myFunction = lambda x, *data: relative_humidity_from_temps(data[0], x) - data[1]
    zGuess = -100
    data = (dry_temp_F, relative_humid)
    out = fsolve(myFunction, zGuess, args=data)

    # out, _ = iter_solver(
    #     -100, relative_humid,
    #     lambda x: relative_humidity_from_temps(
    #         dry_temp_F, x, p_atm_psia=p_atm_psia),
    #     initial_increment=0.001,
    #     num_iters_max=num_iters_max, precision=precision)

    # Wet bulb temperature can't be greater than dry bulb temp:
    return out[0]


def wet_bulb_temperature_empiric(
        dry_temp_F: float, relative_humid: float) -> float:
    """Empiric calculation of the wet bulb temperature for P_atm.

    Ref to (eq1) [http://journals.ametsoc.org/doi/pdf/10.1175/JAMC-D-11-0143.1]
    """
    dry_temp_c = (dry_temp_F - 32) * 5 / 9
    rel_humid_p = relative_humid * 100
    if -2.33 * dry_temp_c + 28.33 > rel_humid_p:  # From Fig 3. Tw Error > 1 ºC
        print('WARNING: The empiric formulation for wetbulb temperature at '
              'this point ({}, {}) is out of range. Expect a '
              'considerable error.'.format(dry_temp_c, relative_humid))

    return (dry_temp_c * atan(0.151977 * sqrt(rel_humid_p + 8.313659))
            + atan(dry_temp_c + rel_humid_p)
            - atan(rel_humid_p - 1.676331)
            + 0.00391838 * (rel_humid_p ** 1.5) * atan(0.023101 * rel_humid_p)
            - 4.686035) * 9 / 5 + 32


# def degree_of_saturation(w_kg_kga, wsat_kg_kga):
#     """Ratio of air humidity ratio to humidity ratio of saturated moist air.
#
#     µ (no dimension), eq (12) 2009 ASHRAE Handbook—Fundamentals (SI).
#     """
#     return w_kg_kga / wsat_kg_kga
#
#
# def relative_humidity(p_vapor_kpa, p_sat_vapor_kpa):
#     """Ratio of the mole fraction of water vapor x_w in a given moist
#     air sample to the mole fraction xws in an air sample
#     saturated at the same temperature and pressure.
#
#     Eq (24) 2009 ASHRAE Handbook—Fundamentals (SI).
#     """
#     return p_vapor_kpa / p_sat_vapor_kpa
