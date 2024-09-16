import Elements as El
from numpy import real_if_close, array
from RealParameters import P0, Pp0, Q0
from lmfit import Parameters, minimize, fit_report
from NumericFunctions import ModSquared
FittingValues = Parameters()

FittingValues.add_many(("KsiAS",5.76311793,True,1,10),
                ("KsiAP",3.65783562,True,1,10),
                ("KsiCS",4.37252536,True,1,10),
                ("KsiCP",3.30091699,True,1,10),
                ("CAS",-0.88133599,True,-1.4,1.4),
                ("CAXv",1.13135549,True,-1.4,1.4),
                ("CCS",1.11176923,True,-1.4,1.4),
                ("CCXv",-0.85612875,True,-1.4,1.4),
                ("CAXc",-0.89456498,True,-1.4,1.4),
                ("CCXc",-1.10740896,True,-1.4,1.4),
                ("CAZv",1.14448738,True,-1.4,1.4),
                ("CCZv",-0.84674614,True,-1.4,1.4),
                )

RealParametersValue ={
    "P0" : El.GaAs().P0,
    "Pp0" : El.GaAs().Pp0,
    "Q0" : El.GaAs().Q0
}
AnionQuant = El.GaAs().AnionPrincipalQuantumNuber
CationQuant = El.GaAs().CationPrincipalQunatumNumber
TargetData = array([RealParametersValue["P0"],real_if_close(RealParametersValue["Pp0"]*(-1j)),RealParametersValue["Q0"],2,2,2,2,0,0,0,0,0])


def Residue(Pars, data):
    ParsValues = Pars.valuesdict()
    
    P0ToFit = P0(Pars, PrincipallQuantumNumberOfAnion=AnionQuant, PrincipalQuantumNumberOfCation=CationQuant).CalculateParameterValue()
    Pp0ToFit = real_if_close(Pp0(Pars, PrincipallQuantumNumberOfAnion=AnionQuant, PrincipalQuantumNumberOfCation=CationQuant).CalculateParameterValue()*(-1j))
    Q0ToFit = Q0(Pars, PrincipallQuantumNumberOfAnion=AnionQuant, PrincipalQuantumNumberOfCation=CationQuant).CalculateParameterValue()
    
    NormalizationOne = ModSquared(ParsValues["CAS"])*ModSquared(ParsValues["CCS"])
    NormalizationTwo = ModSquared(ParsValues["CAXv"])*ModSquared(ParsValues["CCXv"])
    NormalizationThree = ModSquared(ParsValues["CAXc"])*ModSquared(ParsValues["CCXc"])
    NormalizationFour = ModSquared(ParsValues["CAZv"])*ModSquared(ParsValues["CCZv"])
    
    CurrentFittingValue = array([P0ToFit,Pp0ToFit,Q0ToFit,NormalizationOne,NormalizationTwo,NormalizationThree,NormalizationFour,0,0,0,0,0])
    return CurrentFittingValue - data

out = minimize(Residue, FittingValues, kws={'data': TargetData})

print(fit_report(out))