"""Code for particle class"""
from constants import c

class ParticleError(Exception):
    pass

class Particle:
    """Particle class detailing particle's type, position, motion, etc."""
    def __init__(self,particleType,**kwargs):
        self.identifyType(particleType)
        self.position = [0.,0.,0.] #m
        self.momentum = [0.,0.,0.] #MeV
        for key,val in kwargs.items():
            if "pos" in key and len(val)==3:
                self.position = [float(x) for x in val]
            if "mom" in key and len(val)==3:
                self.momentum = [float(p) for p in val]
            if "vel" in key and len(val)==3:
                self.momentum = [v*self.mass/c for v in val]

    def identifyType(self,particleType):
        """Sets particle attributes based on type name"""
        if "neutrino" in particleType or "nu" in particleType:
            if "E" in particleType or "el" in particleType:
                self.type = "nuE"
                self.mass = 0 #MeV
                self.charge = 0 #e
            elif "Mu" in particleType or "mu" in particleType:
                self.type = "nuMu"
                self.mass = 0 #MeV
                self.charge = 0 #e
            else:
                raise ParticleError("Unrecognized particle type "+particleType)
            if "bar" in particleType.lower() or "anti" in particleType:
                self.type += "Bar"
        elif "pi" in particleType:
            if "+" in particleType or "positive" in particleType:
                self.type = "pi+"
                self.mass = 139.57018 #MeV
                self.charge = 1 #e
            elif "-" in particleType or "negative" in particleType:
                self.type = "pi-"
                self.mass = 139.57018 #MeV
                self.charge = -1 #e
            elif "0" in particleType or "neutral" in particleType:
                self.type = "pi0"
                self.mass = 134.9766 #MeV
                self.charge = 0 #e
            else:
                raise ParticleError("Unrecognized particle type "+particleType)
        elif "mu" in particleType:
            if "+" in particleType or "positive" in particleType:
                self.type = "mu+"
                self.mass = 105.658369 #MeV
                self.charge = 1 #e
            elif "-" in particleType or "negative" in particleType:
                self.type = "mu-"
                self.mass = 105.658369 #MeV
                self.charge = -1 #e
            else:
                raise ParticleError("Unrecognized particle type "+particleType)
        elif particleType=="positron" or particleType=="e+":
            self.type = "e+"
            self.mass = 0.5109989 #MeV
            self.charge = 1 #e
        elif particleType=="electron" or particleType=="e-" or particleType=="e":
            self.type = "e-"
            self.mass = 0.5109989 #MeV
            self.charge = -1 #e
        elif particleType=="proton" or particleType=="p+" or particleType=="p":
            self.type = "p+"
            self.mass = 938.27203 #MeV
            self.charge = 1 #e
        elif particleType=="neutron" or particleType=="n0" or particleType=="n":
            self.type = "n0"
            self.mass = 939.56536 #MeV
            self.charge = 0
        else:
            raise ParticleError("Unrecognized particle type "+str(particleType))

