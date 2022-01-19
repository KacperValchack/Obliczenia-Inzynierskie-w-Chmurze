import math

local_longitude_deg = 22.0945515
local_longitude_h = 1.473

#method taken from Orbital Mechanics for Engineering Students
def get_lst_deg (day, month, year, time_local_s):
    j0 = 367*year - math.floor((7*(year+math.floor((month+9)/12)))/4) + math.floor(275*month/9) + day + 1721013.5
    t0 = (j0-2451545)/36525
    theta_g0 = 100.4606184 + 36000.77004*t0 + 0.000387933*t0**2 - 2.583*10**(-8)*t0**3
    if theta_g0 > 360:
        while theta_g0 > 360:
            theta_g0 = theta_g0 - 360
    if theta_g0 < 0:
        while theta_g0 < 0:
            theta_g0 = theta_g0 + 360
    time_greenwich_s = time_local_s - 3600
    time_greenwich_h = time_greenwich_s / 3600
    theta_g = theta_g0 + 360.98564724 * time_greenwich_h/24
    theta = theta_g + local_longitude_deg
    return theta


