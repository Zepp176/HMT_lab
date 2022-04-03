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
C_r = C_c/C_h

T_hi = np.array([54.65, 47.20, 40.30, 42.27, 43.38, 37.23, 40.30, 51.27, 41.57])
T_ci = np.array([18.66, 18.18, 17.50, 17.05, 17.14, 19.97, 17.50, 17.81, 17.12])
T_ho = np.array([37.27, 35.95, 30.55, 31.54, 34.29, 30.60, 26.06, 30.38, 24.19])
T_co = np.array([33.11, 28.55, 25.37, 24.73, 31.39, 28.17, 22.39, 29.12, 26.61])

Q_c = C_c*(T_co - T_ci)
Q_h = C_h*(T_hi - T_ho)
Q = (Q_c + Q_h) / 2

def DTm_co(T_hi, T_ho, T_ci, T_co):
    DT1 = T_hi - T_ci
    DT2 = T_ho - T_co
    return (DT2-DT1)/np.log(DT2/DT1)

def DTm_counter(T_hi, T_ho, T_ci, T_co):
    DT1 = T_ho - T_ci
    DT2 = T_hi - T_co
    return (DT2-DT1)/np.log(DT2/DT1)

U = []
ref = [0, 0, 0, 1, 1, 0, 0, 0, 1]

for i in range(9):
    if ref[i] == 0:
        U.append( Q[i] / ( A[i] * DTm_co(T_hi[i],T_ho[i],T_ci[i],T_co[i]) ) )
    elif ref[i] == 1:
        U.append( Q[i] / ( A[i] * DTm_counter(T_hi[i],T_ho[i],T_ci[i],T_co[i]) ) )

C_min = C_c
NTU = U*A/C_min

def f_co(NTU, C_r):
    return ( 1 - np.exp(-NTU*(1+C_r)) ) / (1+C_r)

def f_counter(NTU, C_r):
    if C_r == 1:
        return NTU/(1+NTU)
    return ( 1 - np.exp(-NTU*(1-C_r)) ) / ( 1 - C_r*np.exp(-NTU*(1-C_r)) )

Q_max = C_min*(T_hi - T_ci)
eps = Q/Q_max

print("Case   m_c   m_h    Q        NTU       U        A       eps")
for i in range(9):
    print("  {}  | {:.0f}  | {:.0f}  | ".format(i+1, m_c[i]*3600, m_h[i]*3600), end="")
    print("{:.2f} | ".format(Q[i]), end="")
    print("{:.4f} | ".format(NTU[i]), end="")
    print("{:07.3f} | {:.3f} | ".format(U[i], A[i]), end="")
    print("{:.4f}".format(eps[i]))

NTU_ref = np.linspace(0, 1.75, 100)
eps_ref = (1 - np.exp(-NTU_ref*2)) / 2
eps2_ref = NTU_ref/(1 + NTU_ref)
eps3_ref = (1 - np.exp(-NTU_ref*(1 - 20/35))) / (1 - 20/35*np.exp(-NTU_ref*(1 - 20/35)) )

eps1 = np.sort(np.concatenate((eps[:3], eps[5:8])))
NTU1 = np.sort(np.concatenate((NTU[:3], NTU[5:8])))
eps2 = np.sort(np.array([eps[3], eps[-1]]))
NTU2 = np.sort(np.array([NTU[3], NTU[-1]]))

plt.figure(figsize=(7, 5))
plt.plot(NTU_ref, eps_ref)
plt.plot(NTU_ref, eps2_ref)
plt.plot(NTU_ref, eps3_ref)
plt.plot(NTU1, eps1, 'ko--')
plt.plot(NTU2, eps2, 'ko--', label='_nolegend_')
plt.plot([NTU[4]], [eps[4]], 'ok', label='_nolegend_')
plt.grid()
plt.title("$\epsilon$ - NTU curves")
plt.xlabel("NTU")
plt.ylabel("$\epsilon$")
plt.ylim(bottom=0, top=0.65)
plt.xlim(left=0)
plt.legend(["co-flow & $C_r=1$", "counter-flow & $C_r=1$", "counter-flow & $C_r=4/7$", "experimental data"])
plt.tight_layout()
plt.savefig("fig1.png", dpi=300)











# Partie mystique d'Antoine


NTU_raw = np.arange(0,1.3,0.01)
plt.figure(figsize=(7, 5))
plt.plot(NTU_raw,f_co(NTU_raw,1),label = 'Theoretical curve at $C_r = 1$')
plt.plot(NTU[0:3],eps[0:3],'o',label = "Empirical values from the lab")
plt.grid()
plt.title(" $\epsilon - NTU$ curve", fontsize = 15)
plt.xlabel("NTU",fontsize = 15)
plt.ylabel('$\epsilon$',fontsize = 15)
plt.legend(fontsize = 15)
plt.show()
  
A_x = 2*np.pi*0.006*np.arange(0,2.1,0.01)
T_h_co = (T_hi[1] + T_ci[1])/2 + 0.5*(T_hi[1] - T_ci[1])*np.exp(-U[1]*(2/C_h[1])*A_x)
T_c_co = (T_hi[1] + T_ci[1])/2 - 0.5*(T_hi[1] - T_ci[1])*np.exp(-U[1]*(2/C_h[1])*A_x)
T_h_counter = np.linspace(T_ho[3],T_hi[3],len(A_x))
T_c_counter = np.linspace(T_ci[3],T_co[3],len(A_x))

plt.figure(figsize=(8, 5))
plt.plot(np.arange(0,2.1,0.01),T_h_co,'r',label = "$T_h$ - coflow case")
plt.plot(np.arange(0,2.1,0.01),T_c_co,'b',label = "$T_c$ - coflow case")
plt.plot(np.arange(0,2.1,0.01),T_h_counter,'r--',label = "$T_h$ - counterflow case")
plt.plot(np.arange(0,2.1,0.01),T_c_counter,'b--',label = "$T_c$ - counterflow case")
plt.xlabel("x [m]",fontsize =15)
plt.ylabel("T [°C]",fontsize =15)
plt.legend(fontsize = 10, loc = 'upper right')
plt.grid()
plt.title("Exponential temperature profiles - coflow and counterflow", fontsize = 15)
plt.show()


dT_co = []
dT_counter = []
for i in range(len(T_c_co)):
    dT_co.append((T_h_co[i]-T_c_co[i]))
    dT_counter.append((T_h_counter[len(T_c_co)-i-1]-T_c_counter[i]))

plt.figure(figsize=(11, 5))
plt.plot(np.arange(0,2.1,0.01)*1/2.1,dT_co,'k--',label = "$\Delta T$ - coflow case")
plt.plot(np.arange(0,2.1,0.01)*1/2.1,dT_counter,'k',label = "$\Delta T$ - counterflow case")
plt.xlabel("x/L",fontsize =15)
plt.ylabel("$\Delta T$ [°C]",fontsize =15)
plt.legend(fontsize = 10, loc = 'upper right')
plt.grid()
plt.title("Profil of the temperature differences between hot/cold flow for coflow and counterflow cases ", fontsize = 15)
plt.show()
