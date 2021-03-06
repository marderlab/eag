# This file defines gating functions as in Liu et al 1998
# also used by O'Leary 2014

# gating functions
boltz(V::Float64,A::Float64,B::Float64) = 1/(1 + exp((V+A)/B))
tauX(V::Float64,A::Float64,B::Float64,D::Float64,E::Float64) = A - B/(1+exp((V+D)/E))

# NaV
mNainf(V::Float64) = boltz(V,25.5,-5.29)
hNainf(V::Float64) = boltz(V,48.9,5.18)
taumNa(V::Float64) = tauX(V,1.32,1.26,120.,-25.)
tauhNa(V::Float64) = (0.67/(1+exp((V+62.9)/-10.0)))*(1.5 + 1/(1+exp((V+34.9)/3.6)))

# CaT
mCaTinf(V::Float64) = boltz(V,27.1,-7.2)
hCaTinf(V::Float64) = boltz(V,32.1,5.5)
taumCaT(V::Float64) = tauX(V,21.7,21.3,68.1,-20.5)
tauhCaT(V::Float64) = tauX(V,105.,89.8,55.,-16.9)

# CaS
mCaSinf(V::Float64) = boltz(V,33.,-8.1)
hCaSinf(V::Float64) = boltz(V,60.,6.2)
taumCaS(V::Float64) = 1.4 + (7/((exp((V+27)/10))+(exp((V+70)/-13))))
tauhCaS(V::Float64) = 60 + (150/((exp((V+55)/9))+(exp((V+65)/-16))))

# A
hAinf(V::Float64) = boltz(V,56.9,4.9)
mAinf(V::Float64) = boltz(V,27.2,-8.7)
taumA(V::Float64) = tauX(V,11.6,10.4,32.9,-15.2)
tauhA(V::Float64) = tauX(V,38.6,29.2,38.9,-26.5)

# KCa
mKCainf(V::Float64,Ca::Float64) = (Ca/(Ca+3))*(1/(1+exp((V+28.3)/-12.6)))
taumKCa(V::Float64) = tauX(V,90.3,75.1,46.,-22.7)

# Kd
mKdinf(V::Float64) = boltz(V,12.3,-11.8)
taumKd(V::Float64) = tauX(V,7.2,6.4,28.3,-19.2)

# H 
mHinf(V::Float64) = boltz(V,70.0,6.0)
taumH(V::Float64) = (272.0 + 1499.0/(1.0+exp((V + 42.2)/-8.73)));

