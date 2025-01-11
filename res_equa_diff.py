import numpy as np


def euler_order_1(DT,T_START, T_MAX, f_p, F0, *args, **kwargs):

    T = [0]

    F =  [F0]
    F_P =  [f_p(F0, *args, **kwargs)]
    i = 0
    while F[-1] > 0 and abs(T[-1]) < abs(T_MAX):

        dt = DT*np.min([1,F[-1]**(3/2)])
        
        F_new = F[-1] + dt*F_P[-1]
        F_P_new= f_p(F[-1])

        F.append(F_new)
        F_P.append(F_P_new)
        T.append(T[-1] + dt)

        i+=1
    F = np.array(F)
    F_P = np.array(F_P)
    T = np.array(T)

    return F[F>0], F_P[F>0], T[F>0]


def RK4(DT,T_START, T_MAX, f_p, F0, *args, **kwargs):

    T = [0]

    F =  [F0]
    F_P =  [f_p(F0, *args, **kwargs)]
    i = 0
    while F[-1] > 0 and abs(T[-1]) < abs(T_MAX):

        k1 = f_p(F[-1], *args, **kwargs)
        k2 = f_p(F[-1] + DT/2*k1, *args, **kwargs)
        k3 = f_p(F[-1] + DT/2*k2, *args, **kwargs)
        k4 = f_p(F[-1] + DT*k3, *args, **kwargs)

        F_new = F[-1] + DT/6*(k1 + 2*k2 + 2*k3 + k4)
        F_P_new= f_p(F[-1])

        F.append(F_new)
        F_P.append(F_P_new)
        T.append(T[-1] + DT)

        i+=1
    F = np.array(F)
    F_P = np.array(F_P)
    T = np.array(T)
    return F[F>0], F_P[F>0], T[F>0]

def cosmo_model(DT, a_pp, A0, A_P0,T_START, T_END, *args, **kwargs):
    T = [T_START]
    a = [A0]
    a_P = [A_P0]
    a_PP = [a_pp(A0, *args, **kwargs)]
    last_a = A0
    i =0

    while last_a > 1e-5 and abs(T[-1]) <T_END :
        dt = DT*np.min([1,a[-1]**(3/2)])
        a_new =  a[-1] + dt*a_P[-1]
        a_P_new = a_P[-1] + dt*a_PP[-1]
        a_PP_new = a_pp(a[-1], *args, **kwargs)

        a.append(a_new)
        a_P.append(a_P_new)
        a_PP.append(a_PP_new)

        T.append( T[-1] + dt)

        last_a = a_new
        i +=1
    return np.array(a),np.array(a_P), np.array(T)






def couple_diff(DT, deltac_pp, deltab_pp, DELTAC_0, DELTAB_0,  DELTAC_P0, DELTAB_P0,T_START, T_END, a, H, rhoc, Omega_m,):

    T = [T_START]
    deltac = [DELTAC_0]
    deltac_P = [DELTAC_P0]
    deltab = [DELTAB_0]
    deltab_P = [DELTAB_P0]
    deltac_PP = [deltac_pp(H[0],a[0], DELTAC_P0,DELTAC_0,DELTAB_0,rhoc[0], Omega_m[0],)]
    deltab_PP = [deltab_pp(H[0],a[0], DELTAB_P0,DELTAB_0,DELTAC_0,rhoc[0], Omega_m[0], )]
    i =0

    while abs(T[-1]) < abs(T_END) and i < len(a):
        dt = DT*np.min([1,a[i]**(3/2)])
        deltac_new =  deltac[-1] + dt*deltac_P[-1]
        deltab_new =  deltab[-1] + dt*deltab_P[-1]

        deltac_P_new = deltac_P[-1] + dt*deltac_PP[-1]
        deltab_P_new = deltab_P[-1] + dt*deltab_PP[-1]

        deltac_PP_new = deltac_pp(H[i], a[i], deltac_P[-1], deltac[-1], deltab[-1],rhoc[i],Omega_m[i], )
        deltab_PP_new = deltab_pp(H[i], a[i], deltab_P[-1], deltab[-1], deltac[-1],rhoc[i],Omega_m[i],)


        deltac.append(deltac_new)
        deltab.append(deltab_new)

        deltac_P.append(deltac_P_new)
        deltab_P.append(deltab_P_new)

        deltac_PP.append(deltac_PP_new)
        deltab_PP.append(deltab_PP_new)

        T.append( T[-1] + dt)
        i +=1

    return np.array(deltac),np.array(deltab)