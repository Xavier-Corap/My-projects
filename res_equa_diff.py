import numpy as np

def euler_order_1(N, T_MAX, f_p, F0):

    T = np.linspace(0,T_MAX,N)
    DT = T_MAX/(N-1)

    F = np.empty_like(T)
    F_P = np.empty_like(T)

    F[0] = F0
    F_P[0] = f_p(F0)
    
    for i in range(N-1):

        F[i+1] = F[i] + DT*f_p(F[i])
        F_P[i+1] = f_p(F[i+1])

    return F, F_P, T



def cosmo_model(DT,a_p, a_pp, A0, T_END, OMEGA_R0=4.6e-5, OMEGA_M0=0.3,OMEGA_L0=0.7, OMEGA_K0=0):
    T = [0]
    a = [A0]
    a_P = [a_p(A0, OMEGA_R0=OMEGA_R0, OMEGA_M0=OMEGA_M0,OMEGA_L0=OMEGA_L0, OMEGA_K0=OMEGA_K0)]
    a_PP = [a_pp(A0, OMEGA_R0=OMEGA_R0, OMEGA_M0=OMEGA_M0,OMEGA_L0=OMEGA_L0, OMEGA_K0=OMEGA_K0)]
    last_a = A0
    i =0

    while last_a > 1e-3 and abs(T[-1]) <T_END :
        dt = DT*np.min([1,a[-1]**(3/2)])
        a_new =  a[-1] + dt*a_P[-1]
        a_P_new = a_P[-1] + dt*a_PP[-1]
        a_PP_new = a_pp(a[-1], OMEGA_R0=OMEGA_R0, OMEGA_M0=OMEGA_M0,OMEGA_L0=OMEGA_L0, OMEGA_K0=OMEGA_K0)

        a.append(a_new)
        a_P.append(a_P_new)
        a_PP.append(a_PP_new)

        T.append( T[-1] + dt)

        last_a = a_new
        i +=1
    return np.array(a), np.array(T)


def euler_order_2(N, T_MAX, f_pp, F0,F_P0, H, a, omega_m, rho_c):
    
    T = np.linspace(0,T_MAX,N)
    DT = T_MAX/(N-1)
    print(DT)
    F = np.empty_like(T)
    F_P = F.copy()
    F_PP = F.copy()

    F[0] = F0
    F_P[0] = F_P0

    for i in range(N-1):

        F_P[i+1] = F_P[i] + DT*f_pp(H[i],a[i],F_P[i], F[i], omega_m[i], rho_c[i])
        F[i+1] = F[i] + DT*F_P[i]
        F_PP[i+1] = f_pp(H[i+1],a[i+1],F_P[i+1], F[i+1], omega_m[i+1], rho_c[i+1])

    return F, F_P, F_PP, T


def RK4(N, T_MAX, f_p, F0, OMEGA_R0=4.6e-5, OMEGA_M0=0.24,OMEGA_L0=0.76, OMEGA_K0=0):

    T = np.linspace(0,T_MAX,N)
    DT = T_MAX/(N-1)

    F = np.empty_like(T)
    F_P = np.empty_like(T)

    F[0] = F0
    F_P[0] = f_p(F0, OMEGA_R0=OMEGA_R0, OMEGA_M0=OMEGA_M0,OMEGA_L0=OMEGA_L0, OMEGA_K0=OMEGA_K0)

    for i in range(N-1):

        k1 = f_p(F[i], OMEGA_R0=OMEGA_R0, OMEGA_M0=OMEGA_M0,OMEGA_L0=OMEGA_L0, OMEGA_K0=OMEGA_K0)
        k2 = f_p(F[i] + DT/2*k1, OMEGA_R0=OMEGA_R0, OMEGA_M0=OMEGA_M0,OMEGA_L0=OMEGA_L0, OMEGA_K0=OMEGA_K0)
        k3 = f_p(F[i] + DT/2*k2, OMEGA_R0=OMEGA_R0, OMEGA_M0=OMEGA_M0,OMEGA_L0=OMEGA_L0, OMEGA_K0=OMEGA_K0)
        k4 = f_p(F[i] + DT*k3, OMEGA_R0=OMEGA_R0, OMEGA_M0=OMEGA_M0,OMEGA_L0=OMEGA_L0, OMEGA_K0=OMEGA_K0)

        F[i+1] = F[i] + DT/6*(k1 + 2*k2 + 2*k3 + k4)
        F_P[i+1] = f_p(F[i+1])

    return F, F_P, T
