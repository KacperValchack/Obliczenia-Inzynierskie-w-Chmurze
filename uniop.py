import math
import find_root
import get_lst_deg


mi = 3.986004419 * 10**14  # gravitational parameter of Earth
Re = 6378000  # Earth equatorial radius in m
flattening = 0.003353  # flattening of the earth
phi_geodetic = math.radians(50.0016071)  # geodetic latitude
h = 325  # altitude


lst_deg = []  # local sideral time in degrees
lst_rad = []  # local sideral time in radians
ra = []  # right ascention in radians
declination = []  # declination in radians
year = []
month = []
day = []
t = []  # time of the observation in seconds
f = open(r"dane_iss_4.txt", "r")


while True:
    line = f.readline()
    words = line.split()
    if not words:
        break
    ra.append(math.radians(float(words[0][0:3]) + float(words[0][4:6])/60 + float(words[0][7:])/3600))
    declination.append(math.radians((-1)*(float(words[1][1:3]) + float(words[1][4:6])/60 + float(words[1][7:])/3600)))
    year.append(float(words[2][6:]))
    month.append(float(words[2][3:5]))
    day.append(float(words[2][0:2]))
    t.append(float(words[3][0:2])*3600 + float(words[3][3:5])*60 + float(words[3][6:]))  # local time in seconds
    lst_deg.append(get_lst_deg.get_lst_deg(day[-1], month[-1], year[-1], t[-1]))
    # print (lst_deg)
    lst_rad.append(lst_deg[-1] * (2*math.pi/360))
    # lst_deg = [44.506, 45.000, 45.499]  # test case 1 parameter! to be commented before calculating preliminary orbit from real data
    # lst_rad = [44.506* (2*math.pi/360), 45.000 * (2*math.pi/360), 45.499 * (2*math.pi/360)]  # test case 1 parameter! to be commented before calculating preliminary orbit from real data
print(f'lst_deg = {lst_deg}\nra = {ra}\ndeclination = {declination}')

f.close()
# lst = (lst_rad[0] + lst_rad[1] + lst_rad[2]) / 3  #mean local sideral time in rad

# observer position vector [in meters]
R = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]
row_nr = 0

while row_nr < 3:
    R[row_nr][0] = (Re/(math.sqrt(1-(2*flattening-flattening*flattening)*math.sin(phi_geodetic)*math.sin(phi_geodetic)))+h)*math.cos(phi_geodetic)*math.cos(lst_rad[row_nr]) #first element of I
    R[row_nr][1] = (Re/(math.sqrt(1-(2*flattening-flattening*flattening)*math.sin(phi_geodetic)*math.sin(phi_geodetic)))+h)*math.cos(phi_geodetic)*math.sin(lst_rad[row_nr])  #first element of J
    R[row_nr][2] = ((Re*(1-flattening)*(1-flattening))/(math.sqrt(1-(2*flattening-flattening*flattening)*math.sin(phi_geodetic)*math.sin(phi_geodetic)))+h)*math.sin(phi_geodetic)  #first element of K
    row_nr += 1

# orbiting body direction cosine vector [dimentionless]
ro = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]
row_nr = 0

while row_nr < 3:
    ro[row_nr][0] = (math.cos(declination[row_nr]) * math.cos(ra[row_nr]))  # first element of I
    ro[row_nr][1] = (math.cos(declination[row_nr]) * math.sin(ra[row_nr]))  # first element of J
    ro[row_nr][2] = (math.sin(declination[row_nr]))  # first element of K
    row_nr += 1
print(f"R = {R}\nro = {ro}")

# time differences calculation [in seconds]
tau = t[2] - t[0]
tau1 = t[0] - t[1]
tau3 = t[2] - t[1]
print(f'tau = {tau}\ntau1 = {tau1}\ntau3 = {tau3}')

# cross product of the observational unit direction
p = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]

p[0][0] = (ro[1][1]*ro[2][2] - ro[1][2]*ro[2][1])  # ro2 x ro3, i component, ro22*ro33-ro23*ro32
p[0][1] = -(ro[1][0]*ro[2][2] - ro[1][2]*ro[2][0])  # ro2 x ro3, j component
p[0][2] = (ro[1][0]*ro[2][1] - ro[1][1]*ro[2][0])  # ro2 x ro3, k component
p[1][0] = (ro[0][1]*ro[2][2] - ro[0][2]*ro[2][1])  # ro1 x ro3, i component
p[1][1] = -(ro[0][0]*ro[2][2] - ro[0][2]*ro[2][0])  # ro1 x ro3, j component
p[1][2] = (ro[0][0]*ro[2][1] - ro[0][1]*ro[2][0])  # ro1 x ro3, k component
p[2][0] = (ro[0][1]*ro[1][2] - ro[0][2]*ro[1][1])  # ro1 x ro2, i component
p[2][1] = -(ro[0][0]*ro[1][2] - ro[0][2]*ro[1][0])  # ro1 x ro2, j component
p[2][2] = (ro[0][0]*ro[1][1] - ro[0][1]*ro[1][0])  # ro1 x ro2, k component
print(f'p = {p}')

# common scalar quantity
D0 = ro[0][0]*p[0][0] + ro[0][1]*p[0][1] + ro[0][2]*p[0][2]
print(f'D0 = {D0}')

# nine scalar quantities
D = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]
r_num = 0
p_num = 0

while r_num < 3:
    while p_num < 3:
        D[r_num][p_num] = R[r_num][0]*p[p_num][0] + R[r_num][1]*p[p_num][1] + R[r_num][2]*p[p_num][2]
        p_num += 1
    p_num = 0
    r_num += 1
print(f'D = {D}')

# scalar position coefficients
A = (1/D0) * ((-1)*D[0][1]*(tau3/tau) + D[1][1] + D[2][1]*(tau1/tau))
B = (1/(6*D0)) * (D[0][1]*(tau3*tau3-tau*tau)*(tau3/tau) + D[2][1]*(tau*tau-tau1*tau1)*(tau1/tau))
E = R[1][0]*ro[1][0] + R[1][1]*ro[1][1] + R[1][2]*ro[1][2]
print(f'A = {A}\nB = {B}\nE = {E}')

# squared scalar distance of the second observation
R2_2 = R[1][0]*R[1][0] + R[1][1]*R[1][1] + R[1][2]*R[1][2]  # ok until here!!!
print(f'R2_2 = {R2_2}')

# coefficients of the scalar distance polynomial for the second observation
a = (-1) * (A*A + 2*A*E + R2_2)
b = (-2) * mi * B * (A + E)
c = (-1) * mi * mi * B * B
print(f"a = {a}, b = {b}, c = {c}")

# the root of the scalar distance polynomial for the second observation
r2 = find_root.find_root_newton(a, b, c)
print(f'r2 = {r2}')  # also ok here (if positive)!!!  around 6800km which is roughly earth radius + iss orbit height

# the slant range
ro1 = (1/D0)*(((6*(D[2][0]*(tau1/tau3)+D[1][0]*(tau/tau3))*r2**3 + mi*D[2][0]*(tau**2-tau1**2)*(tau1/tau3)) / (6*r2**3+mi*(tau**2-tau3**2))) - D[0][0])
ro2 = A + (mi*B)/r2**3
ro3 = (1/D0)*(((6*(D[0][2]*(tau3/tau1)-D[1][2]*(tau/tau1))*r2**3 + mi*D[0][2]*(tau**2-tau3**2)*(tau3/tau1)) / (6*r2**3+mi*(tau**2-tau1**2))) - D[2][2])
print(f'ro1 = {ro1}\nro2 = {ro2}\nro3 = {ro3}')
# print(ro)

# orbiting body position vectors
r = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]
r[0][0] = R[0][0] + ro1 * ro[0][0]  # r1i
r[0][1] = R[0][1] + ro1 * ro[0][1]  # r1j
r[0][2] = R[0][2] + ro1 * ro[0][2]  # r1k
r[1][0] = R[1][0] + ro2 * ro[1][0]  # r2i  position vector!!!
r[1][1] = R[1][1] + ro2 * ro[1][1]  # r2j  position vector!!!
r[1][2] = R[1][2] + ro2 * ro[1][2]  # r2k  position vector!!!
r[2][0] = R[2][0] + ro3 * ro[2][0]  # r3i
r[2][1] = R[2][1] + ro3 * ro[2][1]  # r3j
r[2][2] = R[2][2] + ro3 * ro[2][2]  # r3k
print(f'r = {r}')
# print(f'r2 = {r2}')
# print(f'r2 = {math.sqrt(r[1][0]**2+r[1][1]**2+r[1][2]**2)}')

# Lagrange coefficients
f1 = 1 - 0.5*(mi/r2**3)*tau1**2
f3 = 1 - 0.5*(mi/r2**3)*tau3**2
g1 = tau1 - (1/6)*(mi/r2**3)*tau1**3
g3 = tau3 - (1/6)*(mi/r2**3)*tau3**3
print(f'f1 = {f1}\nf3 = {f3}\ng1 = {g1}\ng3 = {g3}')

# velocity vector for the second observation of the orbiting body
v2 = [0, 0, 0]
v2[0] = (1/(f1*g3-f3*g1))*((-1)*f3*r[0][0]+f1*r[2][0])
v2[1] = (1/(f1*g3-f3*g1))*((-1)*f3*r[0][1]+f1*r[2][1])
v2[2] = (1/(f1*g3-f3*g1))*((-1)*f3*r[0][2]+f1*r[2][2])
print(f'v2 = {v2}')

r1 = math.sqrt(r[0][0]**2+r[0][1]**2+r[0][2]**2)
print(f'r1 = {r1}')
r3 = math.sqrt(r[2][0]**2+r[2][1]**2+r[2][2]**2)
print(f'r2 = {r2}')
print(f'r3 = {r3}')
print(f"Hight of the ISS orbit is {math.sqrt(r[1][0]**2+r[1][1]**2+r[1][2]**2) - Re} meters and it's velocity equals to {math.sqrt(v2[0]**2+v2[1]**2+v2[2]**2)} m/s")

# calculating orbital elements for ISS acording to r2 and v2 vectors
# r[1] = [-6045000, -3490000, 2500000]  # test values only!!! to be removed in the final version!!!
# v2 = [-3457, 6618, 2533]  # test values only!!! to be removed in the final version!!!

# distance
r_scalar = math.sqrt(r[1][0]**2+r[1][1]**2+r[1][2]**2)
print(f'r2 = {r_scalar}')

# velocity
v_scalar = math.sqrt(v2[0]**2+v2[1]**2+v2[2]**2)
print(f'v2 = {v_scalar}')

# radial velocity
vr = (r[1][0]*v2[0] + r[1][1]*v2[1] + r[1][2]*v2[2])/r_scalar
print(f'vr = {vr}')

# specific angular momentum
h = [0, 0, 0]
h[0] = r[1][1]*v2[2] - r[1][2]*v2[1]  # i component
h[1] = (-1) * (r[1][0]*v2[2] - r[1][2]*v2[0])  # j component
h[2] = r[1][0]*v2[1] - r[1][1]*v2[0]  # k component
print(f'h = {h}')

# magnitude of the specific angular momentum
h_scalar = math.sqrt(h[0]**2 + h[1]**2 + h[2]**2)  # specific angular momentum!!!
print(f'h = {h_scalar}')

# inclination
i = math.acos(h[2]/h_scalar)  # inclination !!!
print(f'i = {i}')

# N
K = [0, 0, 1]
N = [0, 0, 0]
N[0] = K[1]*h[2] - K[2]*h[1]  # i component
N[1] = (-1) * (K[0]*h[2] - K[2]*h[0])  # j component
N[2] = K[0]*h[1] - K[1]*h[0]  # k component
print(f'N = {N}')

# magnitude of N
N_scalar = math.sqrt(N[0]**2 + N[1]**2 + N[2]**2)
print(f'N = {N_scalar}')

# RA of the ascending node
if N[1] >= 0:
    omega = math.acos(N[0]/N_scalar)  # RA of the ascending node!!!
else:
    omega = 2*math.pi - math.acos(N[0] / N_scalar)
print(f'capital omega = {omega}')

# eccentricity vector
e = [0, 0, 0]
e[0] = (1/mi)*((v_scalar**2 - mi/r_scalar)*r[1][0] - (r[1][0]*v2[0]+r[1][1]*v2[1]+r[1][2]*v2[2])*v2[0])
e[1] = (1/mi)*((v_scalar**2 - mi/r_scalar)*r[1][1] - (r[1][0]*v2[0]+r[1][1]*v2[1]+r[1][2]*v2[2])*v2[1])
e[2] = (1/mi)*((v_scalar**2 - mi/r_scalar)*r[1][2] - (r[1][0]*v2[0]+r[1][1]*v2[1]+r[1][2]*v2[2])*v2[2])
# e[0] = (1/mi)*((v2[1]*h[2]-v2[2]*h[1])-mi*r[1][0]/r_scalar)
# e[1] = (1/mi)*((v2[2]*h[0]-v2[0]*h[2])-mi*r[1][1]/r_scalar)
# e[2] = (1/mi)*((v2[0]*h[1]-v2[1]*h[0])-mi*r[1][2]/r_scalar)
print(f'e = {e}')

# eccentricity
# e_scalar = 1/mi*(math.sqrt((2*mi-r_scalar*v_scalar*v_scalar)*r_scalar*vr*vr + (mi - r_scalar*v_scalar*v_scalar)**2))  # eccentricity!!!
e_scalar = math.sqrt(e[0]**2 + e[1]**2 + e[2]**2)  # eccentricity!!!
print(f'e = {e_scalar}')

# argument of perigee
if e[2] >= 0:
    omega_small = math.acos((N[0]*e[0] + N[1]*e[1] + N[2]*e[2])/(N_scalar*e_scalar))  # argument of perigee!!!
else:
    omega_small = 2*math.pi - math.acos((N[0] * e[0] + N[1] * e[1] + N[2] * e[2]) / (N_scalar * e_scalar))
print(f'omega small = {omega_small}')

# true anomaly
if vr >= 0:
    theta_scalar = math.acos((e[0]*r[1][0] + e[1]*r[1][1] + e[2]*r[1][2])/(e_scalar*r_scalar))  # true anomaly!!!
else:
    theta_scalar = 2*math.pi - math.acos((e[0]*r[1][0] + e[1]*r[1][1] + e[2]*r[1][2])/(e_scalar*r_scalar))
print(f'theta = {theta_scalar}')

print(f'Orbital elements are as follows:\nh = {h_scalar}\ni = {i}\ncapital omega = {omega}\ne = {e_scalar}\nsmall omega = {omega_small}\ntheta = {theta_scalar}')

#writing results to file:
f = open (r"wyniki.txt", "w")

f.write(f"Hight of the ISS orbit is {math.sqrt(r[1][0]**2+r[1][1]**2+r[1][2]**2) - Re} meters and it's velocity equals to {math.sqrt(v2[0]**2+v2[1]**2+v2[2]**2)} m/s\n"
        f"Orbital elements are as follows:\nh = {h_scalar}\ni = {i}\ncapital omega = {omega}\ne = {e_scalar}\nsmall omega = {omega_small}\ntheta = {theta_scalar}\n")

f.close()
