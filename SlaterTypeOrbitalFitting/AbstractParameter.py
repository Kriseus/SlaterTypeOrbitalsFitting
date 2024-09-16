from numpy import array
from abc import ABC, abstractmethod
from WaveFunctionMethods import SlaterTypeOrbital
from NumericFunctions import ThreeDimIntegral


class Parameter(ABC):
    @abstractmethod
    def __init__(self, Params, PrincipallQuantumNumberOfAnion, PrincipalQuantumNumberOfCation, df = 1e-7, NumOfPts = 5, LatticeConstant = 5.5, hbar = 6.582119569e-16, mass = 0.510998950e6/(2.99792458e18**2)):
        ParamsDict = Params.valuesdict()
        
        self.NumOfPts = NumOfPts
        self.df = df
        self.a=LatticeConstant
        self.hbar=hbar
        self.m0=mass 
        self.RangesOfIntegral = (self.a, self.a, self.a)
        self.PositionCoordination = array([[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]])*self.a/4
        
        self.QuantNumAnion = PrincipallQuantumNumberOfAnion
        self.QuantNumCation = PrincipalQuantumNumberOfCation
        
        
        self.KsiAS=ParamsDict["KsiAS"]
        self.KsiAP=ParamsDict["KsiAP"]
        self.KsiCS=ParamsDict["KsiCS"]
        self.KsiCP=ParamsDict["KsiCP"]
        self.CAS=ParamsDict["CAS"]*1j
        self.CAXv=ParamsDict["CAXv"]
        self.CCS=ParamsDict["CCS"]*1j
        self.CCXv=ParamsDict["CCXv"]
        self.CAXc=ParamsDict["CAXc"]*1j
        self.CCXc=ParamsDict["CCXc"]*1j
        self.CAZv=ParamsDict["CAZv"]
        self.CCZv=ParamsDict["CCZv"]
        
        self.Iteration = 0 
        
    def CalculateAnionSTypeWaveFunction(self, Coordinates):

        x, y, z = Coordinates[0], Coordinates[1], Coordinates[2]

        PsiS = SlaterTypeOrbital(0,0,self.QuantNumAnion,self.KsiAS).CalcWaveFunc(x, y, z)
            
        return PsiS
    def CalculateAnionPTypeWaveFunction(self, Coordinates):
        PsiP = 0
        x, y, z, Direction = Coordinates[0], Coordinates[1], Coordinates[2], Coordinates[3]
        if Direction =="X":
            PsiP += SlaterTypeOrbital(1,1,self.QuantNumAnion,self.KsiAP).CalcWaveFunc(x, y, z)
        else:
            PsiP += SlaterTypeOrbital(1,0,self.QuantNumAnion,self.KsiAP).CalcWaveFunc(x, y, z)
    
        return PsiP
    def CalculateCationSTypeWaveFunction(self, Coordinates):
        x, y, z = Coordinates[0], Coordinates[1], Coordinates[2]
        PsiS = SlaterTypeOrbital(0,0,self.QuantNumCation,self.KsiCS).CalcWaveFunc(x, y, z)
            
        return PsiS
    def CalculateCationPTypeWaveFunction(self, Coordinates):
        PsiP = 0
        x, y, z, Direction = Coordinates[0], Coordinates[1], Coordinates[2], Coordinates[3]
        if Direction =="X":
            PsiP += SlaterTypeOrbital(1,1,self.QuantNumCation,self.KsiCP).CalcWaveFunc(x, y, z)
        else:
            PsiP += SlaterTypeOrbital(1,0,self.QuantNumCation,self.KsiAP).CalcWaveFunc(x, y, z)
            
        return PsiP
    def AnionAnionCoupling(self, x, y, z):
        pass
    def CationCationCoupling(self, x, y, z):
        pass
    def AnionCationCoupling(self, x, y, z):
        
        x1 = x - self.PositionCoordination[self.Iteration][0]
        y1 = y - self.PositionCoordination[self.Iteration][1]
        z1 = z - self.PositionCoordination[self.Iteration][2]
        
        return x1,y1,z1
    def CationAnionCoupling(self, x, y, z):
        
        x1 = x + self.PositionCoordination[self.Iteration][0]
        y1 = y + self.PositionCoordination[self.Iteration][1]
        z1 = z + self.PositionCoordination[self.Iteration][2]
        
        return x1,y1,z1
    def MakeIntegrals(self):
        CationCationIntegral = 0
        AnionAnionIntegral = 0
        AnionCationIntegral = 0 
        CationAnionIntegral = 0
        
        CationCationIntegral += ThreeDimIntegral(self.CationCationCoupling, self.RangesOfIntegral, NumberOfPoints=self.NumOfPts)
        AnionAnionIntegral += ThreeDimIntegral(self.AnionAnionCoupling, self.RangesOfIntegral, NumberOfPoints=self.NumOfPts)
        
        while(self.Iteration<4):
            AnionCationIntegral += ThreeDimIntegral(self.AnionCationCoupling, self.RangesOfIntegral, NumberOfPoints=self.NumOfPts)
            CationAnionIntegral += ThreeDimIntegral(self.CationAnionCoupling, self.RangesOfIntegral, NumberOfPoints=self.NumOfPts)
            self.Iteration += 1
        
        return AnionAnionIntegral, CationCationIntegral, AnionCationIntegral, CationAnionIntegral
    def CalculateParameterValue(self):
        print("zaczynam liczenie!")
        return self.MakeIntegrals()
