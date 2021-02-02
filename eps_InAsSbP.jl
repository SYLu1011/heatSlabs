__precompile__()
module eps_InAsSbP
export G_E1_E2,func_Q,func_Q2,E_T,func_f,InAsSbP_struct, eps_InAsSbP_xy,eps_InAsSbP_imag_xy
#Used this version in the end, accurate imaginary part from the range of 0.1 to 1.5, above 1.5, imaginary part becomes negative"""
#import numpy as np
#import matplotlib.pyplot as plt

e= 1.602*10^(-19)
eV = 1.60218*10^(-19)
eps_0 = 8.8542*10^(-12)
c = 3.0*10^8
hbar = 1.05457*10^(-34) #m2kg/s
hbar_eV = 6.5821*10^(-16) #eV s 
m_e =9.109*10^(-31)#kg
m_0 = 0.9109*10^(-30)
m_0_eV = 0.51099906*10^6 #eV
T = 300.0

#list of constants: list of constants for InAsSbP, AD=InAs, BD=InSb, CD=InP ABC=AsSbP  #Cuevas Paper Table1
#can't find bowing constants for AsSbP
#Bowing Constants for GaIn, AsSb from ZhangThermoConversion
A_InAs = 19.950   #A (eV^1.5) for eps_1, Cuevas
A_InSb = 16.830   #A (eV^1.5) for eps_1, Cueva #A_InP =  6.57     #Need fitting Adachi Table1
C_InAs = 1.750  #C(Dimension-less) for eps_3 Cuevas
C_InSb = 2.64#C_InP = 1.49 #Need fitting Adachi Table1
D_InAs =21.511 #D(Dimension-less) for eps_4 Cuevas
D_InSb = 38.830 #D_InP = 60.4 #Need fitting Adachi Table1
G_E0_InAs = 8.660 #in eV Cuevas
G_E0_InSb = 4.20

@inline function  G_E1_E2(gamma_L,g,T)
    return gamma_L + g*T
end
precompile(G_E1_E2,(Float64,Float64,Float64))
G_E1_InAs = G_E1_E2(0.25,0.3,T)*0.001 #gamma in meV, change to eV #Cuevas Table 3
G_E1_InSb = G_E1_E2(0.15,0.28,T)*0.001
G_E1_InP  = 0.094 # Adachi Table1

G_E2_InAs =  G_E1_E2(403,0.3,T)*0.001 #gamma in meV, change to eV #Cuevas Table 3
G_E2_InSb =  G_E1_E2(633.9,0.49,T)*0.001 
G_E2_InP =  0.10 # Adachi Table1

G_Eg_InAs =  0.350#(eV) Cuevas
G_Eg_InSb = 0.56 
#Cuevas Table2, The energies are in eV, NO bowing constants for E0??Only temperature depence?
E0_InAs = 0.42  #at 300K, E0_T = 1.5,at 0K 1.4
dt0_InAs = 0.28*0.001
b0_InAs = 93.0
E0_InSb = 0.24
dt0_InSb = 0.32*0.001
b0_InSb = 170.0 #E0_InP  = 1.35 #Need fitting Adachi Table1

E0Delta0_InAs = 0.79#-0.42 
dtE0Delta0_InAs = 0.34*0.001
bE0Delta0_InAs  = 248.0
E0Delta0_InSb = 1.20#-0.24
dtE0Delta0_InSb = 0.32*0.001
bE0Delta0_InSb = 170.0 #E0Delta0_GaAs = 1.85#-0.152 dtE0Delta0_GaAs = 0.35*0.001 bE0Delta0_GaAs = 225.0 E0Delta0_GaSb = 1.57#-0.81 dtE0Delta0_GaSb = 0.42*0.001
#E0Delta0_InP = 1.35 + 0.10 # Adachi Table1

E1_InAs = 2.5 #Cuevas Table2 
dt1_InAs = 0.5*0.001
b1_InAs = 0.0
E1_InSb = 2.0 
dt1_InSb = 0.68*0.001
b1_InSb = 132.0
#E1_InP = 3.10 # Adachi Table1
  
E1Delta1_InAs =2.88#-2.61 #Cuevas Table 2
dtE1Delta1_InAs = 0.5*0.001
bE1Delta1_InAs = 0.0
E1Delta1_InSb = 2.49#-2.00
dtE1Delta1_InSb = 0.65*100
bE1Delta1_InSb = 170.0 #E1Delta1_InP = 3.10+0.15 #Adachi  Table1
  
Delta1_InAs = 2.88-2.61 #Cuevas
Delta1_InSb = 2.49-2.00 #Delta1_InP = 0.15 #Adachi  Table1

E2_InAs = 4.70 #from Zhang Thermo Conversion
dt2_InAs = 0.56*0.001
b2_InAs = 0.0
E2_InSb = 3.9
dt2_InSb = 0.54*0.001
b2_InSb = 0.0
E2_InP = 4.7 #Adachi  Table1

Eg_InAs = 1.13 #Cuevas Table 2
dtg_InAs = 0.28*0.001
bg_InAs = 93.0
Eg_InSb = 0.93
dtg_InSb = 0.32*0.001
bg_InSb = 170.0
Eg_InP  = 2.05 #Adachi Table1

#Cuevas paper Table IV in angstrong
alpha_lc_InAs = 6.06
alpha_lc_InSb = 6.48
alpha_lc_InP =  5.87 
#alpha_lc_GaAs = 5.65 alpha_lc_GaSb = 6.10 alpha_lc_GaIn = 0.055  #Cuevas no value #zhang 0.055 alpha_lc_AsSb = -0.023 #Cuevas no value #zhang -0.023

#Fitted Parameters
A_InP =  31.199833249457633
C_InP=  0.5600889435237922
D_InP=  -147.6161185407525
E0_InP = 2.106123911548122
E0Delta0_InP = 1.422040863793526
E1_InP = 3.030856569453735
E1Delta1_InP = 2.88535017475427
Delta1_InP =  -0.3610122449508903
G_E0_InP= 0.12337661235349834
G_Eg_InP=  0.8902957174564856

#Q Parameter in Notes
@inline function  func_Q(x,y,B_AD,B_BD,B_AC,B_CD,C_ABC)
    return x*B_AD + y*B_BD + (1-x-y)*B_CD  + C_ABC*x*y*(1-x-y)
end
precompile(func_Q,(Float64,Float64,Float64,Float64,Float64,Float64,Float64))

 #Q Parameter in Notes, No bowing parameters,  AD=InAs, BD=InSb, CD=InP ABC=AsSbP
@inline function  func_Q2(x,y,B_AD,B_BD,B_CD) 
    return x*B_AD + y*B_BD + (1-x-y)*B_CD 
end
precompile(func_Q2,(Float64,Float64,Float64,Float64,Float64))

@inline function  E_T(E0,dt,beta)
    return E0-dt*300.0^2/(300.0+beta)
end
precompile(E_T,(Float64,Float64,Float64))

@inline function  func_f(z)
    return z^(-2)*(2.0-(1.0+z)^0.5-(1.0-z)^0.5)
end
precompile(func_f,(Float64,))

struct InAsSbPDsc
    x::Float64
    y::Float64
    A_InAsSbP::Float64
    C_InAsSbP::Float64
    D_InAsSbP::Float64
    G_E0_InAsSbP::Float64 #No value found for Gamma_E1
    G_E1_InAsSbP::Float64
    G_Eg_InAsSbP::Float64
    G_E2_InAsSbP::Float64
    E0_InAsSbP::Float64
    E0Delta0_InAsSbP::Float64
    E1_InAsSbP::Float64
    E1Delta1_InAsSbP::Float64
    Delta1_InAsSbP::Float64
    E2_InAsSbP::Float64
    Eg_InAsSbP::Float64
    alpha_lc_InAsSbP::Float64
end

@inline function InAsSbP_struct(x,y)
     #assume equal proportion of As,Sb and P
    #list of constants for InAsSbP, AD=InAs, BD=InSb, CD=InP ABC=AsSbP  #Cuevas Paper Table1
    #Bowing Constants for GaIn, AsSb from ZhangThermoConversion
	E0_T_InAs = E_T(E0_InAs,dt0_InAs,b0_InAs)
	E0_T_InSb = E_T(E0_InSb,dt0_InSb,b0_InSb)
	E0Delta0_T_InAs = E_T(E0Delta0_InAs,dtE0Delta0_InAs,bE0Delta0_InAs)
	E0Delta0_T_InSb = E_T(E0Delta0_InSb,dtE0Delta0_InSb,bE0Delta0_InSb)
	E1_T_InAs = E_T(E1_InAs,dt1_InAs,b1_InAs)
	E1_T_InSb = E_T(E1_InSb,dt1_InSb,b1_InSb)
	E1Delta1_T_InAs = E_T(E1Delta1_InAs,dtE1Delta1_InAs,bE1Delta1_InAs)
	E1Delta1_T_InSb =E_T(E1Delta1_InSb,dtE1Delta1_InSb,bE1Delta1_InSb)        
	E2_T_InAs = E_T(E2_InAs,dt2_InAs,b2_InAs)
	E2_T_InSb = E_T(E2_InSb,dt2_InSb,b2_InSb)
	Eg_T_InAs = E_T(Eg_InAs,dtg_InAs,bg_InAs)
	Eg_T_InSb = E_T(Eg_InSb,dtg_InSb,bg_InSb)
	
	A = func_Q2(x,y,A_InAs,A_InSb,A_InP)
    C= func_Q2(x,y,C_InAs,C_InSb,C_InP)
    D = func_Q2(x,y,D_InAs,D_InSb,D_InP)
    G_E0 = func_Q2(x,y,G_E0_InAs,G_E0_InSb,G_E0_InP) #No value found for Gamma_E1
    G_E1 = func_Q2(x,y,G_E1_InAs,G_E1_InSb,G_E1_InP)
    G_Eg = func_Q2(x,y,G_Eg_InAs,G_Eg_InSb,G_Eg_InP)
    G_E2 = func_Q2(x,y,G_E2_InAs,G_E2_InSb,G_E2_InP)
    E0 = func_Q2(x,y,E0_T_InAs,E0_T_InSb,E0_InP)    #No dt0_InP 
    E0Delta0 = func_Q2(x,y,E0Delta0_T_InAs, E0Delta0_T_InSb,E0Delta0_InP) #no E0Delta0_T_InP
    E1 = func_Q2(x,y,E1_T_InAs,E1_T_InSb,E1_InP) #no E1_T_InP
    E1Delta1 = func_Q2(x,y,E1Delta1_T_InAs,E1Delta1_T_InSb,E1Delta1_InP)
    Delta1= func_Q2(x,y,Delta1_InAs,Delta1_InSb,Delta1_InP) 
    E2= func_Q2(x,y,E2_T_InAs,E2_T_InSb,E2_InP)
    Eg= func_Q2(x,y,Eg_T_InAs,Eg_T_InSb,Eg_InP)
    alpha_lc = func_Q2(x,y,alpha_lc_InAs,alpha_lc_InSb,alpha_lc_InP)
    return InAsSbPDsc(x,y,A,C,D,G_E0,G_E1,G_Eg,G_E2,E0,E0Delta0,E1,E1Delta1, Delta1, E2,Eg ,alpha_lc)
end
precompile(InAsSbP_struct,(Float64,Float64))

@inline function  eps_InAsSbP_xy(E_photon,InAsSbPstruct)
    #GaInAsSb in the energy range 0.5-6eV omega  (hbar*2*math.pi)*c/(0.5*eV) = 2.5um, frequency 
    # m_GaInAsSb = 0.022+0.032*x -0.012*x^2  #Equation 4a) of Zhang thermal Conversion, close to m_0 
    m_InAsSbP = 0.028     #checked online, cannot find 
    eps_1 =InAsSbPstruct.A_InAsSbP*InAsSbPstruct.E0_InAsSbP^(-1.5)*(func_f((E_photon + im*InAsSbPstruct.G_E0_InAsSbP)/InAsSbPstruct.E0_InAsSbP) + 0.5*(InAsSbPstruct.E0_InAsSbP/(InAsSbPstruct.E0Delta0_InAsSbP))^1.5*func_f((E_photon + im*InAsSbPstruct.G_E0_InAsSbP)/InAsSbPstruct.E0Delta0_InAsSbP))
    B1=44*(InAsSbPstruct.E1_InAsSbP+InAsSbPstruct.Delta1_InAsSbP/3)/(InAsSbPstruct.alpha_lc_InAsSbP*InAsSbPstruct.E1_InAsSbP^2)
    B2=44*(InAsSbPstruct.E1_InAsSbP+2*InAsSbPstruct.Delta1_InAsSbP/3)/(InAsSbPstruct.alpha_lc_InAsSbP*(InAsSbPstruct.E1_InAsSbP+InAsSbPstruct.Delta1_InAsSbP)^2)
    eps_2term1 = -B1*((E_photon + im*InAsSbPstruct.G_E1_InAsSbP)/InAsSbPstruct.E1_InAsSbP)^(-2)*log(1-((E_photon + im*InAsSbPstruct.G_E1_InAsSbP)/InAsSbPstruct.E1_InAsSbP)^2)
    eps_2term2 = -B2*((E_photon + im*InAsSbPstruct.G_E1_InAsSbP)/(InAsSbPstruct.E1Delta1_InAsSbP))^(-2)*log(1-((E_photon + im*InAsSbPstruct.G_E1_InAsSbP)/(InAsSbPstruct.E1Delta1_InAsSbP))^2)
    eps_2 = eps_2term1 + eps_2term2 
    eps_3 = InAsSbPstruct.C_InAsSbP*InAsSbPstruct.E2_InAsSbP^2/(InAsSbPstruct.E2_InAsSbP^2-E_photon^2 + im*(-E_photon*InAsSbPstruct.G_E2_InAsSbP))
    #changed from G_E0_GaInAsSb to G_Eg_GaInAsSb
    E_complex = (E_photon + im*InAsSbPstruct.G_Eg_InAsSbP)
    eps_4_1st =-log(InAsSbPstruct.E1_InAsSbP/InAsSbPstruct.Eg_InAsSbP)*InAsSbPstruct.Eg_InAsSbP^2/E_complex^2
    eps_4_2nd = (1+InAsSbPstruct.Eg_InAsSbP/E_complex)^2*log((E_complex+InAsSbPstruct.E1_InAsSbP)/(E_complex+InAsSbPstruct.Eg_InAsSbP))/2
    eps_4_3rd = (1-InAsSbPstruct.Eg_InAsSbP/E_complex)^2*log((E_complex-InAsSbPstruct.E1_InAsSbP)/(E_complex-InAsSbPstruct.Eg_InAsSbP))/2
    eps_4 = (2*InAsSbPstruct.D_InAsSbP/pi)*(eps_4_1st+eps_4_2nd+eps_4_3rd)
    eps = eps_1 + eps_2 + eps_3 + eps_4
    return eps
end
precompile(eps_InAsSbP_xy,(Float64,InAsSbPDsc))

@inline function eps_InAsSbP_imag_xy(enr,InAsSbPstruct)
    #In the code InAs_xSb_yP_1-x-y,  from III-V Ternary and Quaternary Compounds Table 30.14 uses InP_xAs_ySb_1-x-y, and E0 = 0.512+0.030*y-0.183*y^2
    E_bandgap = 0.512 + 0.030*InAsSbPstruct.x-0.183*InAsSbPstruct.x^2   
    if(enr>E_bandgap)
        return imag(eps_InAsSbP_xy(enr,InAsSbPstruct))
    else
        return 0
    end
end
precompile(eps_InAsSbP_imag_xy,(Float64,InAsSbPDsc))
end