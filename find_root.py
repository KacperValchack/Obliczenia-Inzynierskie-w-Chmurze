import math

#newton-raphson method
def find_root_newton(a, b, c):
    r2 = float(7000000)
    while abs((r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c) / (8 * r2 ** 7 + 6 * a * r2 ** 5 + 3 * b * r2 ** 2)) > 0.02:
        # print (r2)
        r2 = r2 - ((r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c) / (8 * r2 ** 7 + 6 * a * r2 ** 5 + 3 * b * r2 ** 2))
        if r2 < 0:  #I don't know if this must be positive. If it must than activate this if statement!!!
            r2 = -r2
        # if abs((r2**8 + a*r2**6 + b*r2**3 + c) / (8*r2**7 + 6*a*r2**5 + 3*b*r2**2)) > 10:
        #     r2 = r2 - ((r2**8 + a*r2**6 + b*r2**3 + c) / (8*r2**7 + 6*a*r2**5 + 3*b*r2**2))
        # elif r2 < 0:
        #     r2 = r2 + 10
        # else:
        #     r2 = r2 - 10
    return r2

#custom method
def find_root_custom(a, b, c):
    changed_direction = False
    dif = [0, 0]
    r2 = float(100000)
    step = 10
    direction = 1
    dif[0] = r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c
    r2 = r2 + direction*step
    dif[1] = r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c
    while abs(r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c) > 0.05:
        #print(f'r2 = {r2}, dif_old = {dif[0]}, dif_new = {dif[1]}')
        if abs(dif[1]) > abs(dif[0]):
            direction = direction * (-1)
            if changed_direction:
                step = step * 0.8
            changed_direction = True
        else:
            changed_direction = False
            step = step * 1.2
        dif[0] = dif[1]
        r2 = r2 + direction * step
        dif[1] = r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c
    return r2



#older unused version of find_root_custom
    # dif[0] = r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c
    # r2 = r2 + 10
    # dif[1] = r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c
    # while abs(r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c) > 0.05:
    #     if abs(dif[1]) > abs(dif[0]):
    #         while abs(dif[1]) > abs(dif[0]):
    #             print(f'r2 = {r2}, dif_old = {dif[0]}, dif_new = {dif[1]}')
    #             r2 = r2 - 10
    #             dif[0] = dif[1]
    #             dif[1] = r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c
    #     else:
    #         while (dif[1]) < abs(dif[0]):
    #             print(f'r2 = {r2}, dif_old = {dif[0]}, dif_new = {dif[1]}')
    #             r2 = r2 + 10
    #             dif[0] = dif[1]
    #             dif[1] = r2 ** 8 + a * r2 ** 6 + b * r2 ** 3 + c
    # return r2
