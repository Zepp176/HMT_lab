import numpy as np
import matplotlib.pyplot as plt

cp = 4184
A1 = 2*np.pi*0.006*2.1
A2 = 2*A1
A = np.array([A1, A1, A1, A1, A1, A1, A2, A2, A2])

m_c = np.array([25, 35, 40, 35, 20, 25, 20, 25, 35])
m_h = np.array([25, 35, 40, 35, 35, 25, 20, 25, 35])

m_c = m_c / 3600
m_h = m_h / 3600

C_c = m_c*cp
C_h = m_h*cp
C_r = C_c/C_h # C_c = C_min & C_h = C_max

T_hi = np.array([54.65, 47.20, 40.30, 42.27, 43.38, 37.23, 40.30, 51.27, 41.57])
T_ci = np.array([18.66, 18.18, 17.50, 17.05, 17.14, 19.97, 17.50, 17.81, 17.12])
T_ho = np.array([37.27, 35.95, 30.55, 31.54, 34.29, 30.60, 26.06, 30.38, 24.19])
T_co = np.array([33.11, 28.55, 25.37, 24.73, 31.39, 28.17, 22.39, 29.12, 26.61])

Q_c = C_c*(T_co - T_ci)
Q_h = C_h*(T_hi - T_ho)

U_c = Q_c/(A*(T_hi - T_ci))
U_h = Q_h/(A*(T_hi - T_ci))

C_min = C_c
NTU_c = U_c*A/C_min
NTU_h = U_h*A/C_min

def f_co(NTU, C_r):
    return ( 1 - np.exp(-NTU*(1+C_r)) ) / (1+C_r)

def f_counter(NTU, C_r):
    if C_r == 1:
        return NTU/(1+NTU)
    return ( 1 - np.exp(-NTU*(1-C_r)) ) / ( 1 - C_r*np.exp(-NTU*(1-C_r)) )

Q_max = C_min*(T_hi - T_ci)

eps_c = np.empty(len(NTU_c))
eps_h = np.empty(len(NTU_h))

for i in range(len(NTU_c)):
    eps_c[i] = f_co(NTU_c[i], C_r[i])
    eps_h[i] = f_co(NTU_h[i], C_r[i])

Q_NTU_c = Q_max*eps_c
Q_NTU_h = Q_max*eps_h

print("Case   m_c   m_h    Q_c      Q_h     NTU_c    NTU_h     U_c       U_h      A       eps_c    eps_h")
for i in range(9):
    print("  {}  | {:.0f}  | {:.0f}  | ".format(i+1, m_c[i]*3600, m_h[i]*3600), end="")
    print("{:.2f} | {:.2f} | ".format(Q_c[i], Q_h[i]), end="")
    print("{:.4f} | {:.4f} | ".format(NTU_c[i], NTU_h[i]), end="")
    print("{:07.3f} | {:06.2f} | {:.3f} | ".format(U_c[i], U_h[i], A[i]), end="")
    print("{:.4f} | {:.4f}".format(eps_c[i], eps_h[i]))

plt.plot(Q_c, "o")
plt.plot(Q_h, "o")
plt.plot(Q_NTU_c, "o")
plt.plot(Q_NTU_h, "o")
plt.ylim(bottom=0)
plt.legend(["cold", "hot", "NTU cold", "NTU hot"])
plt.ylabel("Heat exchanged")
plt.show()



"""
plt.plot(Q_c, "o")
plt.plot(Q_h, "o")
plt.legend(["cold", "hot"])
plt.ylim(bottom=0)
plt.ylabel("Heat exchanged")

plt.show()

plt.plot(U_c, "o")
plt.plot(U_h, "o")
plt.legend(["cold", "hot"])
plt.ylim(bottom=0)
plt.ylabel("U")

plt.show()"""