
# To Calculate Operating Temperature for Flash Distillation of a Two-Component Saturated Liquid Mixture Assuming Real Gas Behaviour. 
# The mixture taken is - Benzene and Ethanol.

# Algorithm For Program

# 1. Guess a value of the operating temperature.
#    a. Guess values of the Rachford Rice Parameter 'K' for each of the components using ideal gas law data.
#    b. Define the Rachford Rice equation based on the guess value of Parameter 'K'. We name this equation as Equation 1. 
#    c. Hence using the above Equation 1 and the guess values of Parameter 'K' taken, obtain the guess value of 'Psi'. The value is termed as guess since its obtained using guess value of 'K'.
#    d. Hence obtain Yi and Xi. These are termed as guess values toosince they are obtained based on the guess values of the Parameter 'K'.
#    e. We know the Guessed operating temperature, guessed Yi and Xi values, guessed Parameter 'K' and Guess value of 'Psi'.Hence we obtain the calculated value of Patrameter 'K' using Gamma-Phi approach.  
#    f. We again define a new Rachford Rice Equation based on the calculated value of the parameter 'K'. Hence we obtain the calculated value of 'Psi'.
#    g. We calculate Residual value of 'Psi'. If value of residual is between [-0.0001,+0.0001] then we stop and jump to step-b.
#    h. If 'Psi' residual value is not in desired bracket range, then we reiterate from step-a using the calculated values of parameter 'K' as our new guess value. We keep doing this till we get 'Psi' in permissible range.
#    i. Once we get residual value of 'Psi' in the permissible range we move to step-.
# 2. We define a new equation 3 which is very similar to Rachford Rice Equation but based on enthalpy balance.
#    a. Using the calculated value of Parameter 'K' and 'Psi', we solve equation 3 to get calculated value of temperature.
#    b. In turn, we calculate the residual value of temperature based on the calculated value and the initial guess value taken.
#       If the temperature residual is in the range [-0.0001,0.0001] then our calculated temperature value is our solution. Hence westop. If not, we move to next step.
#    c. We take calculated temperature value as our new guess and reiterate from step-a till the temperature residual is obtainedin the permissible limit.  

import scipy
from scipy.optimize import fsolve
import numpy as np

# Given Problem Solved Using Gamma-Phi Approach
# Equation of State Used is Redlich Kwong

# Data for Given Sum
F=1 # Total inlet Feed of Saturated Liquid Solution
zA=0.2 # Mole Fraction of Methanol in Inlet Saturated Liquid Feed Solution
zB=1-zA # Mole Fration of Ethanol in Inlet Saturated Liquid Feed Solution
Tin=401.36 # Temperature of Inlet Saturated Liquid Feed Solution in Kelvin
P=150000 # Flash Point Pressure in Pascal
Tref=353.3 #Kelvin , Higher Value of the Boiling Point of Either of the Two Component used in mixture in pure state
R=8.314472

# Antoinne Coefficients for Benzene
aA=13.7819
bA=2726.81
cA=273.16-217.572 

# Antoinne Coefficients for Ethanol
aB=16.8958
bB=3795.17
cB=42.242

# Data for Benzene designated as 'A'

# Defining Three Functions to Calculate Saturation Pressure of Pure Vapour, Specific Heat Capacity of Pure Vapour and Specific Heat Capacity of Pure Liquid at a given Temperature for Pure Benzene

# Saturation Pressure of Benzene at a Given Temperature  
def PsatA(T):
    Psat=1000*scipy.exp((aA)-((bA)/(T-cA))) # Where T is in Kelvin and PsatA is in Pascal
    return Psat

# Specific Heat Capacity of Benzene at a Given Temperature in Ideal Gas State 
def CpvA(T):
    Cp=R*(-0.206+(0.039064*T)-(13.301*(T**2)*(10**-6))) # Where T is in Kelvin and CpvA is in Joule/mol.Kelvin
    return Cp

# Specific Heat Capacity of Benzene at a Given Temperature in Liquid State
def CplA(T):
    Cp=R*(-0.747+(0.06796*T)-(37.78*(10**-6)*(T**2))) # Where T is in Kelvin and CplA is in Joule/mol.Kelvin
    return Cp

# Critical Properties of Benzene
PcA=48.98*10**5 # in Pascal
TcA=562.2 # Kelvin

# Data for Ethanol designated as 'B'

# Defining Three Functions to Calculate Saturation Pressure of Pure Vapour, Specific Heat Capacity of Pure Vapour and Specific Heat Capacity of Pure Liquid at a given Temperature for Pure Ethanol

# Saturation Pressure of Ethanol at a Given Temperature
def PsatB(T):
    Psat=1000*scipy.exp((aB)-((bB)/(T-cB))) # Where T is in Kelvin and PsatB is in Pascal
    return Psat
    
# Specific Heat Capacity of Ethanol at a Given Temperature in Ideal Gas State 
def CpvB(T):
    Cp=R*(3.518+(0.020001*T)-(6.002*(T**2)*(10**-6))) # Where T is in Kelvin and CpvB is in Joule/mol.Kelvin
    return Cp
    
# Specific Heat Capacity of Ethanol at a Given Temperature in Liquid State
def CplB(T):
    Cp=R*(33.866-(0.1726*T)+(349.17*(10**-6)*(T**2))) # Where T is in Kelvin and CplB is in Joule/mol.Kelvin
    return Cp
    
# Critical Properties of Ethanol
PcB=61.48*10**5 # Pascal
TcB=513.9 # Kelvin


# Defining a Matrix Function of Benzene and Ethanol, for Integral of Specific Heat Capacity in Vapour Phase with Temperature for limits - Tref to T 
# This Function Provides Us With Corresponding Integral Value of the respective Gas in Real State 
def intCpV(T):
    f1=(R*((-0.206*(T-Tref))+(0.5*0.039064*(T**2-Tref**2))-((1/3)*13.301*(T**3-Tref**3)*(10**-6))))+(3*a(Pc,Tc)[0]*scipy.log(Volg[0]/(Volg[0]+b(Pc,Tc)[0]))*((T**-0.5)-(Tref**-0.5))/(2*b(Pc,Tc)[0]))
    f2=(R*((3.518*(T-Tref))+(0.5*0.020001*(T**2-Tref**2))-((1/3)*6.002*(T**3-Tref**3)*(10**-6))))+(3*a(Pc,Tc)[1]*scipy.log(Volg[1]/(Volg[1]+b(Pc,Tc)[1]))*((T**-0.5)-(Tref**-0.5))/(2*b(Pc,Tc)[1]))
    f=scipy.array([f1,f2])
    return f

# Defining a Matrix Function of Benzene and Ethanol, for Integral of Specific Heat Capacity in Liquid Phase with Temperature for limits - Tref to T 
def intCpL(T):
    f1=R*((-0.747*(T-Tref))+(0.5*0.06796*(T**2-Tref**2))-((1/3)*37.78*(10**-6)*(T**3-Tref**3)))
    f2=R*((33.866*(T-Tref))-(0.5*0.1726*(T**2-Tref**2))+((1/3)*349.17*(10**-6)*(T**3-Tref**3)))
    f=scipy.array([f1,f2])
    return f

# Creating an Array of Critical Temperature and Critical Pressure Values Respectively for Benzene-Ethanol System    
Pc=scipy.array([PcA,PcB])
Tc=scipy.array([TcA,TcB])    

# EOS Parameter Calculation by Forming a Matrix
def a(Pc,Tc):
    a1=(0.42748*R**2*Tc[0]**2.5)/Pc[0]
    a2=(0.42748*R**2*Tc[1]**2.5)/Pc[1]
    a=scipy.array([a1,a2])
    return a
def b(Pc,Tc):
    b1=(0.08664*R*Tc[0])/Pc[0]
    b2=(0.08664*R*Tc[1])/Pc[1]
    b=scipy.array([b1,b2])
    return b

# Calculations

# Calculation of Dew Point for Given Sytem of Benzene-Ethanol
def TdpCalc(w):
    x=w[0]
    T=w[1]
    f1=(zA*P)-(x*PsatA(T))
    f2=(zB*P)-((1-x)*PsatB(T))
    f=scipy.array([f1,f2])
    return f
wGuess=scipy.array([0.1,330])
w=fsolve(TdpCalc,wGuess)
Tdp=w[1]
print('The Dew Point Temperature for Inlet Feed in Kelvin is ')
print(Tdp)

# Calculation of Bubble Point for Given System of Benzene-Ethanol
def TbpCalc(v):
    y=v[0]
    T=v[1]
    f1=(y*P)-(zA*PsatA(T))
    f2=((1-y)*P)-((zB)*PsatB(T))
    f=scipy.array([f1,f2])
    return f
vGuess=scipy.array([0.1,330])
v=fsolve(TbpCalc,vGuess)
Tbp=v[1]
print('The Bubble Point Temperature for Inlet Feed in Kelvin is ')
print(Tbp)

# We know that our Operating Temperature Lies Between the Dew Point and Bubble Point

Tg=Tbp + 0.0001 # Guess Value for Operating Temperature
yAg=0.747
yBg=1-yAg
xAg=0.916
xBg=1-xAg
kAg=yAg/xAg # For Benzene, we obtain guess value of Parameter 'K' of Rachford Rice Equation
kBg=(1-yAg)/(1-xAg) # For Ethanol, we obtain guess value of Parameter 'K' of Rachford Rice Equation
Tres=0.2 # A residual value between Guess Temperature and Calculated Temperature is Assumed Initially

while (Tres>0.0001 or Tres<0.0001):
    kAg=kAg # For Benzene
    kBg=kBg # For Ethanol
    
    # Here We Define The Rachford Rice Equation Function based on Guess Value of Parameter 'K'
    def eq1(psi):
        eq1=(zA*(kAg-1)/((psi*kAg)+1-psi))+(zB*(kBg-1)/((psi*kBg)+1-psi))
        return eq1
    psiGuess=0.9
    psig=fsolve(eq1,psiGuess)
    Vg=F*psig
    Lg=F-Vg
    xAg=(F*zA)/(Vg*kAg+Lg)
    yAg=xAg*kAg
    xBg=1-xAg
    yBg=1-yAg
    
    # RK EOS Defined for Saturation Pressure    
    def RKsat(Vol):
        f1=PsatA(Tg)-((R*Tg)/(Vol[0]-b(Pc,Tc)[0]))+(a(Pc,Tc)[0]/((Tg**0.5)*Vol[0]*(Vol[0]+b(Pc,Tc)[0])))
        f2=PsatB(Tg)-((R*Tg)/(Vol[1]-b(Pc,Tc)[1]))+(a(Pc,Tc)[1]/((Tg**0.5)*Vol[1]*(Vol[1]+b(Pc,Tc)[1])))
        f=scipy.array([f1,f2])
        return f
    
    # RK EOS Defined for Operating Pressure
    def RK(Vol):
        f1=P-((R*Tg)/(Vol[0]-b(Pc,Tc)[0]))+(a(Pc,Tc)[0]/((Tg**0.5)*Vol[0]*(Vol[0]+b(Pc,Tc)[0])))
        f2=P-((R*Tg)/(Vol[1]-b(Pc,Tc)[1]))+(a(Pc,Tc)[1]/((Tg**0.5)*Vol[1]*(Vol[1]+b(Pc,Tc)[1])))
        f=scipy.array([f1,f2])
        return f

    VolGuesssat=scipy.array([0.9999*R*Tg/PsatA(Tg),0.9999*R*Tg/PsatB(Tg)])
    Volgsat=fsolve(RKsat,VolGuesssat)

    VolGuess=scipy.array([0.9999*R*Tg/P,0.9999*R*Tg/P])
    Volg=fsolve(RK,VolGuess)
    
    # Function To Calculate Fugacity Coefficient of Pure Vapour Phase under Saturated Conditions
    def PurePHIsat(Vol): # Calculated at Saturation Pressure
        f1=scipy.exp((PsatA(Tg)*Vol[0]/(R*Tg))-1-scipy.log((PsatA(Tg)*Vol[0]/(R*Tg)))+scipy.log(Vol[0]/(Vol[0]-b(Pc,Tc)[0]))+(a(Pc,Tc)[0]*scipy.log(Vol[0]/(Vol[0]+b(Pc,Tc)[0]))/(b(Pc,Tc)[0]*R*Tg**1.5)))
        f2=scipy.exp((PsatB(Tg)*Vol[1]/(R*Tg))-1-scipy.log((PsatB(Tg)*Vol[1]/(R*Tg)))+scipy.log(Vol[1]/(Vol[1]-b(Pc,Tc)[1]))+(a(Pc,Tc)[1]*scipy.log(Vol[1]/(Vol[1]+b(Pc,Tc)[1]))/(b(Pc,Tc)[1]*R*Tg**1.5)))
        f=scipy.array([f1,f2])
        return f
    
    # Mixing Rules Used As Given in Prausitz
    aMix=((yAg*a(Pc,Tc)[0]**0.5)+(yBg*a(Pc,Tc)[1]**0.5))**2
    bMix=(yAg*b(Pc,Tc)[0])+(yBg*b(Pc,Tc)[1])
    VolMix=(yAg*Volg[0])+(yBg*Volg[1])
    print('The volume based on mixing rules is')
    print(VolMix)
    
    # Function To Calculate Fugacity Coefficient of Component in Mixture Under Operating Pressure Conditions
    def MixPHI(Vol):
        f1=scipy.exp(((b(Pc,Tc)[0]/bMix)*((P*Vol[0]/(R*Tg))-1))-scipy.log(P*(Vol[0]-bMix)/(R*Tg))+((aMix*scipy.log((Vol[0]+bMix*(1+2**0.5))/(Vol[0]+bMix*(1-(2**0.5))))*((b(Pc,Tc)[0]/bMix)-(2*(a(Pc,Tc)[0]/aMix)**0.5))/((2**1.5)*bMix*R*Tg))))
        f2=scipy.exp(((b(Pc,Tc)[1]/bMix)*((P*Vol[1]/(R*Tg))-1))-scipy.log(P*(Vol[1]-bMix)/(R*Tg))+((aMix*scipy.log((Vol[1]+bMix*(1+2**0.5))/(Vol[1]+bMix*(1-(2**0.5))))*((b(Pc,Tc)[1]/bMix)-(2*(a(Pc,Tc)[1]/aMix)**0.5))/((2**1.5)*bMix*R*Tg))))
        f=scipy.array([f1,f2])
        return f
    
    # Calculating Activity Coefficients using UniQuac Theory
    def gamma(xAg):
        rA=3.19;rB=2.11
        qA=2.4;qB=1.97
        qAdash=2.4;qBdash=0.89
        aAB=242.53;aBA=-75.13
        p11=(xAg*rA/(xAg*rA+xBg*rB))
        p12=(xBg*rB/(xAg*rA+xBg*rB))
        p1=scipy.array([p11,p12])
        t11=(xAg*qA/(xAg*qA+xBg*qB))
        t12=(xBg*qB/(xAg*qA+xBg*qB))
        t1=scipy.array([t11,t12])
        t21=(xAg*qAdash/(xAg*qAdash+xBg*qBdash))
        t22=(xBg*qBdash/(xAg*qAdash+xBg*qBdash))
        t2=scipy.array([t21,t22])
        tAB=scipy.exp(-aAB/Tg)
        tBA=scipy.exp(-aBA/Tg)
        t=scipy.array([tAB,tBA])
        z=10
        l11=((z/2)*(rA-qA))-(rA-1)
        l12=((z/2)*(rB-qB))-(rB-1)
        l1=scipy.array([l11,l12])
        gamma1=scipy.exp((scipy.log(p1[0]/xAg))+(z*qA*scipy.log(t1[0]/p1[0])/2)+(p1[1]*(l1[0]-(rA*l1[1]/rB)))-(qAdash*scipy.log(t2[0]+t2[1]*t[1]))+(t2[1]*qAdash*((t[1]/(t2[0]+t2[1]*t[1]))-(t[0]/(t2[1]+t2[0]*t[0])))))
        gamma2=scipy.exp((scipy.log(p1[1]/xBg))+(z*qB*scipy.log(t1[1]/p1[1])/2)+(p1[0]*(l1[1]-(rB*l1[0]/rA)))-(qBdash*scipy.log(t2[1]+t2[0]*t[0]))+(t2[0]*qBdash*((t[0]/(t2[1]+t2[0]*t[0]))-(t[1]/(t2[0]+t2[1]*t[1])))))
        gamma=scipy.array([gamma1,gamma2])
        return gamma

    kAc=gamma(xAg)[0]*PurePHIsat(Volgsat)[0]*PsatA(Tg)*scipy.exp(78.11*(P-PsatA(Tg))/(R*Tg*876000))/(P*MixPHI(Volg)[0])
    kBc=gamma(xAg)[1]*PurePHIsat(Volgsat)[1]*PsatB(Tg)*scipy.exp(46.06844*(P-PsatB(Tg))/(R*Tg*789000))/(P*MixPHI(Volg)[1])

    # Here We Define The Rachford Rice Equation Function based on Guess Value of Parameter 'K'
    def eq2(psi):
        eq2=(zA*(kAc-1)/((psi*kAc)+1-psi))+(zB*(kBc-1)/((psi*kBc)+1-psi))
        return eq2
    psiGuess1=0.5
    psic=fsolve(eq2,psiGuess1)

    Vc=F*psic # Calculated Value of Vapour Molar Flow Rate
    Lc=F-Vc # Calculated Value of Liquid Molar Flow Rate
    
    xAc=(F*zA)/(Vc*kAc+Lc) # Calculated Value of Liquid Mole Fraction of Benzene
    xBc=1-xAc # Calculated Value of Liquid Mole Fraction of Ethanol
    yAc=kAc*xAc # Calculated Value of Vapour Mole Fraction of Benzene
    yBc=1-yAc # Calculated Value of Vpour Mole Fraction of Ethanol
    
    psires=(psig-psic)/psig
    print('psires')
    print(psires)

    if (psires>0.0001 or psires<-0.0001):
        kAg=kAc
        kBg=kBc
        Tg=Tg
    else:
        # An Equation Similar to Rachford Rice Equation is Used involving Entalpy Balance and Excluding Pressure Balance, but Involving the Same Derivation essence as that for Rachford Rice Equation
        def eq3(T):
            f=((zA*(kAc-1)*(intCpL(Tin)[0]))/((psic*kAc*(intCpV(T)[0]))+((1-psic)*(intCpL(Tin)[0]))))+((zB*(kBc-1)*(intCpL(Tin)[1]))/((psic*kBc*(intCpV(T)[1]))+((1-psic)*(intCpL(Tin)[1]))))
            return f
        Tgg=Tg
        Tcalc=fsolve(eq3,Tgg)
        Tres=(Tg-Tcalc)/Tg
        print('Tres')
        print(Tres)
        if (Tres>0.0001 or Tres<-0.0001):
            Tg=Tg + 0.0001
            kAg=kAg
            kBg=kBg
            
        else:
            print('The operating temperature is')
            print(Tg)
print('Hello')
            
    





