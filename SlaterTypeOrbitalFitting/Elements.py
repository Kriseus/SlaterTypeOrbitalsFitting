from abc import ABC, abstractmethod

class ChemicalBound(ABC):
    @abstractmethod
    def __init__(self):
        self.CationPrincipalQunatumNumber = None 
        self.AnionPrincipalQuantumNuber = None 
        self.P0 = None
        self.Pp0 = None 
        self.Q0 = None

class GaAs(ChemicalBound):
    def __init__(self):
        super().__init__()
        self.CationPrincipalQunatumNumber = 4
        self.AnionPrincipalQuantumNuber = 4
        self.P0 = 9.343
        self.Pp0 = -0.509
        self.Q0 = 8.350

        
class GaN(ChemicalBound):
    def __init__(self):
        super().__init__()
        self.CationPrincipalQunatumNumber = 4
        self.AnionPrincipalQuantumNuber = 2
        self.P0 = 7.511
        self.Pp0 = 2.597
        self.Q0 = 10.308

        
class InSb(ChemicalBound):
    def __init__(self):
        super().__init__()
        self.CationPrincipalQunatumNumber = 5
        self.AnionPrincipalQuantumNuber = 5
        self.P0 = 8.553
        self.Pp0 = -0.414 
        self.Q0 = 7.573
