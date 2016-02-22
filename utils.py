from math import sqrt, pi, log10, log, exp
from array import array

class DummyClass():
    pass

def deltaPhi(phi1, phi2):
  PHI = abs(phi1-phi2)
  if PHI<=pi:
      return PHI
  else:
      return 2*pi-PHI

def deltaR(eta1, phi1, eta2, phi2):
  deta = eta1-eta2
  dphi = deltaPhi(phi1,phi2)
  return sqrt(deta*deta + dphi*dphi)

def matching(eta,phi,offJet_eta,offJet_phi,offJet_num):
    dRMin = 99.
    index = -1
    for i in range(0,offJet_num):
        dR = deltaR(eta,phi,offJet_eta[i],offJet_phi[i])
        if dR < 0.5 and dR < dRMin:
            dRMin = dR
            index = i
    
    return index
    

def SetVariable(tree,name,option='F',lenght=1,maxLenght=100):
    if option is 'F': arraytype='f'
    elif option is 'O': arraytype='i'
    elif option is 'I': arraytype='i'
    else:
        print 'option ',option,' not recognized.'
        return

    if not type(lenght) is str:
        maxLenght = lenght
        lenght = str(lenght)
    variable = array(arraytype,[0]*maxLenght)
    if maxLenght>1: name = name + '['+lenght+']'
    tree.Branch(name,variable,name+'/'+option)
    return variable

