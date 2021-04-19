# Python Program for the manuscript "True Nulls-Free Wireless Power Transfer 
# Based on Magnetoinductive Waveguide with Alternate Coupling Polarities"
# submitted to IEEE Trans. on Power Electronics.
#
# Program: SingleZm_LinearACMIW, Version 1.0, 18 April 2021
# Linear ACMI line with M cells
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
L = 2.647e-6
R = 1.055
w0 = 2*cmath.pi*13.56e6
Q = w0*L/R
k = 0.09108                       # Nominal Rx coupling at 65-mm height
kRx55 = 0.119768
kRx65 = k
kRx75 = 0.070197
Zo = w0*k*L                       # Characteristic impedance at w0 and no loss
X = 1j*Zo                         # Positive Mutual impedance between adjacent cells
alphaD = cmath.asinh(1/(2*k*Q))
betaD = cmath.pi/2                # Assumming operation at resonance
GammaD = alphaD + 1j*betaD
RL = Zo                           # Set Rx's load resistance to Zo 
#Zs = Zo                           # Initially set Zs to Zo, then changed to average Zin
Zs = 19.0820                      # average Zin: use this for recalculation efficiency 
                                  # and normalized power
Zm = (w0*kRx65*L)**2/(RL+R)       # Change kRx- for different height here
ZLeff = Zm*RL/(RL+R)              # Effective power delivered to RL

# Settings of coupling polarties and Zmi in vector form
# Zcell = [Zm1 Zm2 ... ZM] <-- modify fixed impedance at each cell here
#                              for reactive termination and multi-loads
# Note: Termination cell loaded with +/-X for ACMIW, X=0 for typical MIW
Zcell = [0,0,0,0,0,0,-X]

# ucell = [u1 u2 ... uM] <-- modify coupling polarities here
ucell = [-1,+1,-1,+1,-1,+1,+1]    # Coupling polarities for cells 1 to N
                                  # @ZM = -X, ucellM is always +1
M = len(Zcell)

# Performance calculations for Zm at cell 1 to cell M
#   - Definition of Parameters
#   subscript "cw" for clockwise variables
#   subscript "ccw" for counterclockwise variables
#   rhoS: Reflection at the source
#   rho: Reflection
#   tau: transmittance
#   T: total transmittance upto cell n
#   Icell_norm: Normalized loop current at each cell
#   Zin,inn: Input impedance
#   Zcell_loaded: Zcell with load
#   Eff: Power transfer efficiency
#   PLav: Load power normalized to available source power

rhoS = 0
rho = [1]*M
tau = [1]*(M+1)
T = [1]*(M+1)
Icell_norm = [0]*M
Icellnorm = [0]*M
ZeT1 = 0
Zin = 0
Zinn = [0]*M
Zcell_loaded = Zcell
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
# T: Transmittance at cell m
    Icellnorm = T*cmath.exp(-m*GammaD)*(1+rho)/(1+rhoS)
    return(Icellnorm)


def Zeff(rho,X,GammaD):
# Effective termination impedance calculation
    ZeT = X*(cmath.exp(-GammaD) + rho*cmath.exp(GammaD))/(1 + rho)
    return(ZeT)


def Linear(Zcellloaded,uX,X,GammaD):
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


def powercalc(rhoS,Icellnorm,X,GammaD,Zs,R,ZLeff,q):
# Calculating Powers
    VSource = 1
    Pav = VSource**2/(4*Zs)
    ZeT1 = Zeff(rhoS,X,GammaD)
    Zin = ZeT1 + R
    ISource = VSource/(Zs + Zin)
    PSource = Zin.real*abs(ISource)**2
    PLoad = ZLeff.real*abs(ISource*Icellnorm[q])**2
    Eff = PLoad/PSource
    PLav = PLoad/Pav
    return(Eff, PLav, Zin)

################################################################################
# Main program
print("Performance of Linear ACMI line with", M, "cells")
print("@ Parameter values")
print("Cell's self-inductance =", L)
print("Cell's total loss      =", R)
print("Coupling factor        =", k)
print("Chacteristic impedance =", Zo)
print("---------------------------------------------------- ")

for q in range (0,M,1):
    Zcell_loaded[q] = Zcell[q] + Zm
    outputLinear = Linear(Zcell_loaded,ucell,X,GammaD)
    rhoS = outputLinear[0]
    Icell_norm = outputLinear[1]
    outputpowercalc = powercalc(rhoS,Icell_norm,X,GammaD,Zs,R,ZLeff,q)
    Eff[q] = outputpowercalc[0] 
    PLav[q] = outputpowercalc[1]
    Zinn[q] = outputpowercalc[2]
    print("ZL at cell",q+1)
    print(" Efficiency =",Eff[q])
    print(" Normalized Load Power =", PLav[q])
    print("---------------------------------------------------- ")
    Zcell_loaded[q] = 0

# Put this value for Zs and recalculation after first run
Zinavg = sum(Zinn)/M
print("Average Zin =", Zinavg)

