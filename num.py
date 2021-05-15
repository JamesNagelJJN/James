import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import linalg as sla
import pandas as pd
import openpyxl
csfont = {'fontname': 'Times New Roman'}
D = 0
neigs = 9

A0 = 16
V0 = 74.9
a0 = 0.6
r0 = 0.97
R0 = r0*A0**(1/3)

A1 = 40
V1 = 51.4
a1 = 0.76
r1 = 1.28
R1 = r1*A1**(1/3)
xmin = -10
xmax = 30
n = 200
x = np.linspace(xmin, xmax, n)
d = R1*3.5 - D
T = []
delta = x[1] - x[0]
print(delta, d)


def oxy(x):
    return -V0/(1+np.exp((abs(x)-R0)/a0))


def cal(x):
    return -V1/(1+np.exp((abs(x-d)-R1)/a1))


def coob():
    cob = []
    for i in x:
        cob.append(oxy(i) + cal(i))
    return cob


def schrodinger1d(n):
    H = sparse.eye(n,n, format='lil')*-2
    K = coob()
    for i in range(n-1):
        if i <= n-2:
            H[i, i+1] = 1
        if i >= 0:
            H[i, i-1] = 1
        H[i, i] += (delta**2) * (K[i]) * 0.04818696
    (evl, eig) = sla.eigs(H, k=neigs, which='SR')
    for i in range(neigs):
        eig[:, i] = eig[:, i] / np.sqrt(np.trapz(np.conj(eig[:, i])*eig[:, i], x))
    return eig, evl/0.04818696


pp = 51

Q = []
www = [R1*3, R1*2, R1*1.5, 0]
for i in range(4):
    d = www[i]
    Q.append(coob())

EE = {'0': Q[0],
      '1': Q[1],
      '2': Q[2],
      '3': Q[3],
      '4': x
      }
df = pd.DataFrame(EE, columns=['0','1','2','3','4'])
df.to_excel(r'C:\Users\James\Desktop\plottttttttttttttt.xlsx', index = False, header=True)

fig, axs = plt.subplots(4)
for i in range(4):
    axs[i].plot(x, Q[i])
    plt.text(10, -200, 'x (fm)', ha='center')
    plt.text(-19, 190, 'V (MeV)', va='center', rotation='vertical')
plt.show()

Q = []
d = R1*3.5 - D

plt.ion()
fig, axs = plt.subplots(2)
Y = []
for i in range(pp):
    SE = schrodinger1d(n)
    Y = SE[0].real
    Q.append(SE[1].real)
    D += 0.3
    d = R1*3.5 - D
    T.append(d)
    plt.xlabel('x (fm)')
    plt.ylabel('$V_{TCSM}$ (MeV)')
    for j in range(neigs):
        O = []
        for i in range(n):
            O.append(Y[i][j] * Y[i][j])
        for i in range(n):
            O[i] += j
        axs[1].plot(x,O)
    axs[0].plot(x, coob())
    plt.draw()
    plt.pause(0.00001)
    axs[0].clear()
    axs[1].clear()
plt.clf()
plt.show()
plt.ioff()
plt.rcParams.update({'font.size': 22})
P = []
print(Q)
plt.plot(T,Q)
plt.show()
PP = []
for j in range(neigs):
    WW = []
    for i in range(len(T)):
        WW.append(Q[i][j])
    PP.append(WW)
print(PP[0])
EE = {'0': T,
      '1': PP[0],
      '2': PP[1],
      '3': PP[2],
      '4': PP[3],
      '5': PP[4],
      '6': PP[5],
      '7': PP[6],
      '8': PP[7],
      '9': PP[8],
      }

df = pd.DataFrame(EE, columns=['0','1','2','3','4','5','6','7','8','9'])
df.to_excel(r'C:\Users\James\Desktop\plotttttt.xlsx', index = False, header=True)
