#parameters 
k1, km1, k2, km2 
Kcat, KM 
k3, k4, k5 

#species = [L, R, LR, LRa, Tfn, TFc, ptfc]
dL/dt    = -k1·L·R + k₋1·LR + k₋2·LRa
dR/dt    = -k1·L·R + k₋1·LR + k₋2·LRa
dLR/dt   =  k1·L·R - k₋1·LR - k2·LR
dLRa/dt  =  k2·LR - k₋2·LRa
dTFc/dt  = -kcat·LRa·TFc/(KM + TFc) + k5·TFcn
dpTFc/dt =  kcat·LRa·TFc/(KM + TFc) - k3·pTFcn
dpTFcn/dt =  k3·pTFc - k4·TFcn
dTFcn/dt  =  k4·pTFcn - k5·TFcn