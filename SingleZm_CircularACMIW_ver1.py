# Python Program for the manuscript "True Nulls-Free Wireless Power Transfer 
# Based on Magnetoinductive Waveguide with Alternate Coupling Polarities"
# submitted to IEEE Trans. on Power Electronics.
#
# Program: SingleZm_LinearACMIW, Version 1.0, 18 April 2021
# 1D-Linear ACMI line with M cells
# Single loaded cell with the load Zm (=Zo) moving from cell 1 to cell M
# By Apisak Worapishet
#    Chaowalit Rukluea
#    Sarawuth Chaimool
#    Yan Zhao
# Copyright (R), 2021.

import math
import cmath

# ACMI Line's electrical parameters
# See Table VI in the manuscript for parameters definitions
L = 2.440e-6
R = 0.986
w0 = 2*cmath.pi*13.56e6
Q = w0*L/R
k = 0.09267                     # Nominal Rx coupling at 65-mm height
kRx55 = 0.119768
kRx65 = k
kRx75 = 0.070197
Zo = w0*k*L                     # Characteristic impedance at w0 and no loss
X = 1j*Zo                       # Positive Mutual impedance between adjacent cells
alphaD = cmath.asinh(1/(2*k*Q))
betaD = cmath.pi/2              # Assumming operation at resonance
GammaD = alphaD + 1j*betaD
RL = Zo                         # Set Rx's load resistance to Zo
#Zs = Zo + 1i*2*Zo;              # Initially set Zs to conjugate if Zo, then changed 
                                # to conjugate of average Zin
Zs = 22.4574 + 1j*35.5804       # Conjugate of average Zin: use this for 
                                # recalculation efficiency and normalized power
Zm = (w0*kRx65*L)**2/(RL+R)     # Change kRx- for different height here
ZLeff = Zm*RL/(RL+R)            # Effective power delivered to RL

# Settings of coupling polarties and Zmi in vector form
# Zcell = [Zm1 Zm2 ... ZM] <-- modify fixed impedance at each cell here
#                              for reactive termination and multi-loads
# ucell = [u1 u2 ... uM] <-- modify coupling polarities here
Zcell_cw = [0,0,0,0,0,0,0]
Zcell_ccw = [0,0,0,0,0,0,0]
ucell = [-1,+1,-1,+1,-1,+1,-1]  # Coupling polarities for cells 1 to N
M = len(Zcell_cw)

# Performance calculations for Zm at cell 1 to cell M
#   - Definition of Parameters
#   subscript "cw" for clockwise variables
#   subscript "ccw" for counterclockwise variables
#   rhoS: Reflection at the source
#   rho: Reflection
#   tau: transmittance
#   T: total transmittance upto cell n
#   Icell_norm: Normalized loop current at each cell
#   ZeT1: Effective terminating impedance at cell 1
#   ZTr1: Transfer impedance from cell 1 to cell M
#   Zin: Input impedance
#   Zcell_loaded: Zcell with load
#   Eff: Power transfer efficiency
#   PLav: Load power normalized to available source power

Zcell_loaded_cw = Zcell_cw
Zcell_loaded_ccw = Zcell_ccw
rhoS_cw = 0
rhoS_ccw = 0
rho = [1]*M
rho_cw = [1]*M
rho_ccw = [1]*M
tau = [1]*(M+1)
tau_cw = [1]*(M+1)
tau_ccw = [1]*(M+1)
T = [1]*(M+1)
ZeT1 = [0]*M
ZeT1_cw = [0]*M
ZeT1_ccw = [0]*M
ZTr1_cw = [0]*M
ZTr1_ccw = [0]*M
Zcw1 = [0]*M
Zccw1 = [0]*M
Zin = [0]*M
Itemp = [0]*M
Icell_norm = 0
Icellnorm_total = [0]*M
Icellnorm_cw = [0]*M
Icellnorm_ccw = [0]*M
ISource = [0]*M
PSource = [0]*M
PLoad = [0]*M
Eff = [0]*M
PLav = [0]*M


def TauRho(Zm,rhop,X,GammaD):
    NUM=((1+rhop*cmath.exp(-2*GammaD))*Zm + rhop*cmath.exp(-2*GammaD)*X*(cmath.exp(GammaD) - cmath.exp(-GammaD)))
    DEM=((1+rhop*cmath.exp(-2*GammaD))*Zm - X*(cmath.exp(GammaD) - cmath.exp(-GammaD)))
    rho = -NUM/DEM
    tau = (1 + rho)/(1 + rhop*cmath.exp(-2*GammaD))
    return(tau, rho)

def round_complex(x):
    return complex(round(x.real,4),round(x.imag,4))

def Itot(m,T,rho,rhoS,GammaD):
# Normalized total loop current at cell m
    Icell_norm = T*cmath.exp(-m*GammaD)*(1+rho)/(1+rhoS)
    return(Icell_norm)

def Zeff(rho,X,GammaD):
# Effective terminating impedance calculation
    ZeT = X*(cmath.exp(-GammaD) + rho*cmath.exp(GammaD))/(1 + rho)
    return(ZeT)

def Linear(Icellnorm,Zcellloaded,uX,X,GammaD):
# rhoS: Reflection at the source
# Initialization
    N = len(Zcellloaded)   # Total number of cells
    rhop = cmath.exp(+2*GammaD)
    tau[0] = 1
# Calculating rho and tau at cell m
    for m in range(N-1,-1,-1):
        outputTauRho = TauRho(Zcellloaded[m],rhop,X,GammaD)
        tau[m+1] = outputTauRho[0] 
        rho[m] = outputTauRho[1]
        rhop = rho[m]
# Calculating total current at cell m (wrt Is = I1 = +1)
    rhoS = rho[0]
    T[0] = tau[0]
    for m in range (0,N,1):     # Is = I1 = +1 is assumed
        Icellnorm[m] = Itot(m,T[m],rho[m],rhoS,GammaD)
        T[m+1] = T[m]*uX[m]*tau[m+1]
    return(rhoS, Icellnorm)


################################################################################
# Main program
print("Performance of Circular ACMI line with", M, "cells")
print("@ Parameter values")
print("Cell's self-inductance =", L)
print("Cell's total loss      =", R)
print("Coupling factor        =", k)
print("Chacteristic impedance =", Zo)
print("---------------------------------------------------- ")

for q in range (0,M,1):
    Zcell_loaded_cw[q] = Zcell_cw[q] + Zm;
    outputLinear_cw = Linear(Icellnorm_cw,Zcell_loaded_cw,ucell,X,GammaD)
    rhoS_cw = outputLinear_cw[0]

    if q == 0:
        outputLinear_ccw = Linear(Icellnorm_ccw,Zcell_loaded_ccw,ucell,X,GammaD)
        rhoS_ccw = outputLinear_ccw[0]
        Icellnorm_ccw = outputLinear_ccw[1]
    else:
        Zcell_loaded_ccw[M-q] = Zcell_ccw[q]+Zm;
        outputLinear_ccw = Linear(Icellnorm_ccw,Zcell_loaded_ccw,ucell,X,GammaD)
        rhoS_ccw = outputLinear_ccw[0]
 
# Calculating effctive terminating and trans- impedances at cell 1
    ZeT1_cw[q] = Zeff(rhoS_cw,X,GammaD)
    ZTr1_cw[q] = ucell[M-1]*X*Icellnorm_cw[M-1]
    Zcw1[q] = ZeT1_cw[q] + ZTr1_cw[q]
    ZeT1_ccw[q] = Zeff(rhoS_ccw,X,GammaD)
    ZTr1_ccw[q] = ucell[M-1]*X*Icellnorm_ccw[M-1]
    Zccw1[q] = ZeT1_ccw[q] + ZTr1_ccw[q]
    ZeT1[q] = Zcw1[q] + Zccw1[q]
    
# Calculating cw and ccw currents then summing for total currents
    Itemp = Icellnorm_ccw[::-1]
    Icellnorm_ccw[1:] = Itemp[:M-1]
    for m in range (0,M,1):
        Icellnorm_total[q] = Icellnorm_ccw[q] + Icellnorm_cw[q]
    Icellnorm_total[0] = 1                   # normalized totoal current to cell 1

# Calculating efficiency and normalized power
    VSource = 1
    Pav = VSource**2/(4*Zs.real)
    Zin[q] = R + ZeT1[q]
    ISource[q] = VSource/(Zs + Zin[q])
    PSource[q] = Zin[q].real*abs(ISource[q])**2
    PLoad[q] = ZLeff.real*abs(ISource[q]*Icellnorm_total[q])**2
    Eff[q] = PLoad[q]/PSource[q]
    PLav[q] = PLoad[q]/Pav
    print("ZL at cell",q+1)
    print(" Efficiency =",Eff[q])
    print(" Normalized Load Power =", PLav[q])
    print("----------------------------------------------------")

# Reset load impedance for next round of calculation
    Zcell_loaded_cw[q] = 0
    if q > 0:
        Zcell_loaded_ccw[M-q] = 0;

# Put this value for Zs and recalculation after first run
Zinavg = sum(Zin)/M
print("Average Zin =", Zinavg)

