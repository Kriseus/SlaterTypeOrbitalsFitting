
from numpy import conjugate, real_if_close
from NumericFunctions import Derivative
from AbstractParameter import Parameter

class P0(Parameter):
    def __init__(self, Params, PrincipallQuantumNumberOfAnion=4, PrincipalQuantumNumberOfCation=4, df=1e-7, NumOfPts=10):
        super().__init__(Params, PrincipallQuantumNumberOfAnion, PrincipalQuantumNumberOfCation, df, NumOfPts)
        
    def CalculateAnionSTypeWaveFunction(self, Coordinates):
        return super().CalculateAnionSTypeWaveFunction(Coordinates)         
    def CalculateAnionPTypeWaveFunction(self, Coordinates):
        return super().CalculateAnionPTypeWaveFunction(Coordinates)   
    def CalculateCationSTypeWaveFunction(self, Coordinates):
        return super().CalculateCationSTypeWaveFunction(Coordinates)
    def CalculateCationPTypeWaveFunction(self, Coordinates):
        return super().CalculateCationPTypeWaveFunction(Coordinates)
    
    def AnionAnionCoupling(self, x, y, z):
        super().AnionAnionCoupling(x, y, z)
        return conjugate(self.CalculateAnionSTypeWaveFunction([x,y,z]))*Derivative(self.CalculateAnionPTypeWaveFunction, [x,y,z,"X"], "x", self.df)
    def CationCationCoupling(self, x, y, z):
        super().CationCationCoupling(x, y, z)
        return conjugate(self.CalculateCationSTypeWaveFunction([x,y,z]))*Derivative(self.CalculateCationPTypeWaveFunction, [x,y,z,"X"], "x", self.df)
    def AnionCationCoupling(self, x, y, z):
        x1, y1, z1 =super().AnionCationCoupling(x, y, z)

        return conjugate(self.CalculateAnionSTypeWaveFunction([x,y,z]))*Derivative(self.CalculateCationPTypeWaveFunction, [x1,y1,z1,"X"],"x",self.df)
    def CationAnionCoupling(self, x, y, z):
        x1, y1, z1=super().CationAnionCoupling(x, y, z)
        return conjugate(self.CalculateCationSTypeWaveFunction([x,y,z]))*Derivative(self.CalculateAnionPTypeWaveFunction, [x1,y1,z1,"X"],"x",self.df)
    
    def MakeIntegrals(self):
        return super().MakeIntegrals()
    
    def CalculateParameterValue(self):
        AnionAnion, CationCation, AnionCation, CationAnion = super().CalculateParameterValue()
        
        AnionAnion = 0.5 * conjugate(self.CAS) * self.CAXv * (-1j * self.hbar) * AnionAnion
        CationCation = 0.5 * conjugate(self.CCS) * self.CCXv * (-1j * self.hbar) * CationCation
        AnionCation = 0.5 * conjugate(self.CAS) * self.CCXv * (-1j * self.hbar) * AnionCation
        CationAnion = 0.5 * conjugate(self.CCS) * self.CAXv * (-1j * self.hbar) * CationAnion
        print("skonczylem liczenie")
        return real_if_close((AnionAnion + CationCation + AnionCation + CationAnion)*(self.hbar/self.m0))

    
class Pp0(Parameter):
    def __init__(self, Params, PrincipallQuantumNumberOfAnion=4, PrincipalQuantumNumberOfCation=4, df=1e-7, NumOfPts=10):
        super().__init__(Params, PrincipallQuantumNumberOfAnion, PrincipalQuantumNumberOfCation, df, NumOfPts)
        
    def CalculateAnionSTypeWaveFunction(self, Coordinates):
        return super().CalculateAnionSTypeWaveFunction(Coordinates)          
    def CalculateAnionPTypeWaveFunction(self, Coordinates):
        return super().CalculateAnionPTypeWaveFunction(Coordinates)   
    def CalculateCationSTypeWaveFunction(self, Coordinates):
        return super().CalculateCationSTypeWaveFunction(Coordinates)  
    def CalculateCationPTypeWaveFunction(self, Coordinates):
        return super().CalculateCationPTypeWaveFunction(Coordinates)
    
    def AnionAnionCoupling(self, x, y, z):
        super().AnionAnionCoupling(x, y, z)
        return conjugate(self.CalculateAnionSTypeWaveFunction([x,y,z]))*Derivative(self.CalculateAnionPTypeWaveFunction, [x,y,z,"X"], "x", self.df)
    def CationCationCoupling(self, x, y, z):
        super().CationCationCoupling(x, y, z)
        return conjugate(self.CalculateCationSTypeWaveFunction([x,y,z]))*Derivative(self.CalculateCationPTypeWaveFunction, [x,y,z,"X"], "x", self.df)
    def AnionCationCoupling(self, x, y, z):
        x1, y1, z1 = super().AnionCationCoupling(x, y, z)
        
        return conjugate(self.CalculateAnionSTypeWaveFunction([x,y,z]))*Derivative(self.CalculateCationPTypeWaveFunction, [x1,y1,z1,"X"],"x",self.df)
    def CationAnionCoupling(self, x, y, z):
        x1, y1, z1 = super().CationAnionCoupling(x, y, z)
        
        return conjugate(self.CalculateCationSTypeWaveFunction([x,y,z]))*Derivative(self.CalculateAnionPTypeWaveFunction, [x1,y1,z1,"X"],"x",self.df)  

    def MakeIntegrals(self):
        return super().MakeIntegrals()
    
    def CalculateParameterValue(self):
        AnionAnion, CationCation, AnionCation, CationAnion = super().CalculateParameterValue()
        
        AnionAnion = 0.5 * conjugate(self.CAS) * self.CAXv * (-1j * self.hbar) * AnionAnion
        CationCation = 0.5 * conjugate(self.CCS) * self.CCXv * (-1j * self.hbar) * CationCation
        AnionCation = 0.5 * conjugate(self.CAS) * self.CCXv * (-1j * self.hbar) * AnionCation
        CationAnion = 0.5 * conjugate(self.CCS) * self.CAXv * (-1j * self.hbar) * CationAnion
        print("skonczylem liczenie")
        return real_if_close((AnionAnion + CationCation + AnionCation + CationAnion)*(self.hbar/self.m0))

class Q0(Parameter):
    def __init__(self, Params, PrincipallQuantumNumberOfAnion=4, PrincipalQuantumNumberOfCation=4, df=1e-7, NumOfPts=10):
        super().__init__(Params, PrincipallQuantumNumberOfAnion, PrincipalQuantumNumberOfCation, df, NumOfPts)
                
    def CalculateAnionPTypeWaveFunction(self, Coordinates):
        return super().CalculateAnionPTypeWaveFunction(Coordinates)   
    def CalculateCationPTypeWaveFunction(self, Coordinates):
        return super().CalculateCationPTypeWaveFunction(Coordinates)      
    
    def AnionAnionCoupling(self, x, y, z):
        super().AnionAnionCoupling(x, y, z)
        return conjugate(self.CalculateAnionPTypeWaveFunction([x,y,z,"Z"]))*Derivative(self.CalculateAnionPTypeWaveFunction, [x,y,z,"X"], "y", self.df)
    def CationCationCoupling(self, x, y, z):
        super().CationCationCoupling(x, y, z)
        return conjugate(self.CalculateCationPTypeWaveFunction([x,y,z,"Z"]))*Derivative(self.CalculateCationPTypeWaveFunction, [x,y,z,"X"], "y", self.df)
    def AnionCationCoupling(self, x, y, z):
        x1, y1, z1 = super().AnionCationCoupling(x, y, z)
        
        return conjugate(self.CalculateAnionPTypeWaveFunction([x,y,z,"Z"]))*Derivative(self.CalculateCationPTypeWaveFunction, [x1,y1,z1,"X"],"y",self.df)
    def CationAnionCoupling(self, x, y, z):
        x1, y1, z1 = super().CationAnionCoupling(x, y, z)

        return conjugate(self.CalculateCationSTypeWaveFunction([x,y,z,"Z"]))*Derivative(self.CalculateAnionPTypeWaveFunction, [x1,y1,z1,"X"],"y",self.df)      
    def MakeIntegrals(self):
        return super().MakeIntegrals()
    def CalculateParameterValue(self):
        AnionAnion, CationCation, AnionCation, CationAnion = super().CalculateParameterValue()
        
        AnionAnion = 0.5 * conjugate(self.CAXc) * self.CAZv * (-1j * self.hbar) * AnionAnion
        CationCation = 0.5 * conjugate(self.CCXc) * self.CCZv * (-1j * self.hbar) * CationCation
        AnionCation = 0.5 * conjugate(self.CAXc) * self.CCZv * (-1j * self.hbar) * AnionCation
        CationAnion = 0.5 * conjugate(self.CCXc) * self.CAZv * (-1j * self.hbar) * CationAnion
        print("skonczylem liczenie")
        return real_if_close((AnionAnion + CationCation + AnionCation + CationAnion)*(self.hbar/self.m0))
    
    