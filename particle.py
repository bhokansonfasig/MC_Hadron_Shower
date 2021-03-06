"""Code for particle class"""
from itertools import count
from numpy import sqrt,sin,cos,arccos,arctan
from constants import c, pi
idCounter = count()

class ParticleError(Exception):
    pass

class Particle:
    """Particle class detailing particle's type, position, motion, etc."""
    def __init__(self,particleType,**kwargs):
        self.identifyType(particleType)
        self.id = next(idCounter)
        self.position = [0.,0.,0.] #m
        self.momentum = [0.,0.,0.] #MeV
        for key,val in kwargs.items():
            if key=="id":
                self.id = val
            elif "pos" in key and len(val)==3:
                self.position = [float(x) for x in val]
            elif "mom" in key:
                if isinstance(val,list) and len(val)==3:
                    self.momentum = [float(p) for p in val]
                elif isinstance(val,float) or isinstance(val,int):
                    self.buildMomentumFromScalar(val,"momentum",kwargs)
            elif "vel" in key:
                if isinstance(val,list) and len(val)==3:
                    self.momentum = [v*self.mass/c for v in val]
                elif isinstance(val,float) or isinstance(val,int):
                    self.buildMomentumFromScalar(val,"velocity",kwargs)
            elif key=="energy" or key=="E":
                self.buildMomentumFromScalar(val,"energy",kwargs)
            elif "kinetic" in key.lower() or key.lower()=="ke":
                self.buildMomentumFromScalar(val,"kinetic",kwargs)
            else:
                if key!="theta" and key!="phi":
                    print("Warning: argument",key,"in Particle not recognized")


    def identifyType(self,particleType):
        """Sets particle attributes based on type name"""
        if "neutrino" in particleType or "nu" in particleType:
            if "E" in particleType or "el" in particleType:
                self.type = "nuE"
                self.mass = 0 #MeV
                self.charge = 0 #e
                self.lifetime = None
            elif "Mu" in particleType or "mu" in particleType:
                self.type = "nuMu"
                self.mass = 0 #MeV
                self.charge = 0 #e
                self.lifetime = None
            elif "Tau" in particleType or "tau" in particleType:
                self.type = "nuTau"
                self.mass = 0 #MeV
                self.charge = 0 #e
                self.lifetime = None
            else:
                raise ParticleError("Unrecognized particle type "+particleType)
            if "bar" in particleType.lower() or "anti" in particleType:
                self.type += "Bar"
        elif "pi" in particleType:
            if "+" in particleType or "positive" in particleType:
                self.type = "pi+"
                self.mass = 139.57018 #MeV
                self.charge = 1 #e
                self.lifetime = 2.6033e-8 #s
            elif "-" in particleType or "negative" in particleType:
                self.type = "pi-"
                self.mass = 139.57018 #MeV
                self.charge = -1 #e
                self.lifetime = 2.6033e-8 #s
            elif "0" in particleType or "neutral" in particleType:
                self.type = "pi0"
                self.mass = 134.9766 #MeV
                self.charge = 0 #e
                self.lifetime = 8.4e-17 #s
            else:
                raise ParticleError("Unrecognized particle type "+particleType)
        elif "mu" in particleType:
            if "+" in particleType or "positive" in particleType:
                self.type = "mu+"
                self.mass = 105.658369 #MeV
                self.charge = 1 #e
                self.lifetime = 2.1969811e-6 #s
            elif "-" in particleType or "negative" in particleType:
                self.type = "mu-"
                self.mass = 105.658369 #MeV
                self.charge = -1 #e
                self.lifetime = 2.1969811e-6 #s
            else:
                raise ParticleError("Unrecognized particle type "+particleType)
        elif particleType=="positron" or particleType=="e+":
            self.type = "e+"
            self.mass = 0.5109989 #MeV
            self.charge = 1 #e
            self.lifetime = None
        elif particleType=="electron" or particleType=="e-" or particleType=="e":
            self.type = "e-"
            self.mass = 0.5109989 #MeV
            self.charge = -1 #e
            self.lifetime = None
        elif particleType=="proton" or particleType=="p+" or particleType=="p":
            self.type = "p+"
            self.mass = 938.27203 #MeV
            self.charge = 1 #e
            self.lifetime = None
        elif particleType=="neutron" or particleType=="n0" or particleType=="n":
            self.type = "n0"
            self.mass = 939.56536 #MeV
            self.charge = 0
            self.lifetime = 881.5 #s
        elif "carbon" in particleType.lower() or particleType=="C":
            if "14" in particleType:
                self.type = "C-14"
                self.mass = 13043.94 #MeV
                self.charge = 6
                self.lifetime = None
            else:
                self.type = "C"
                self.mass = 11177.93 #MeV
                self.charge = 6
                self.lifetime = None
        elif "nitrogen" in particleType.lower() or particleType=="N":
            if "16" in particleType:
                self.type = "N-16"
                self.mass = 14909.59 #MeV
                self.charge = 7
                self.lifetime = 10.29 #s
            else:
                self.type = "N"
                self.mass = 13047.2 #MeV
                self.charge = 7
                self.lifetime = None
        elif "oxygen" in particleType.lower() or particleType=="O":
            if "14" in particleType:
                self.type = "O-14"
                self.mass = 13048.92 #MeV
                self.charge = 8
                self.lifetime = 101.85 #s
            elif "15" in particleType:
                self.type = "O-15"
                self.mass = 13975.27 #MeV
                self.charge = 8
                self.lifetime = 176.4 #s
            else:
                self.type = "O"
                self.mass = 14903.3 #MeV
                self.charge = 8
                self.lifetime = None
        elif "fluorine" in particleType.lower() or particleType=="F":
            if "16" in particleType:
                self.type = "F-16"
                self.mass = 14914.59 #MeV
                self.charge = 9
                self.lifetime = 1.674e-19 #s
            else:
                self.type = "F"
                self.mass = 17696.9 #MeV
                self.charge = 9
                self.lifetime = None
        elif "chlorine" in particleType.lower() or particleType=="Cl":
            if "40" in particleType:
                self.type = "Cl-40"
                self.mass = 37232.2 #MeV
                self.charge = 17
                self.lifetime = 117 #s
            else:
                self.type = "Cl-35"
                self.mass = 32573.28 #MeV
                self.charge = 17
                self.lifetime = None
        elif "argon" in particleType.lower() or particleType=="Ar":
            self.type = "Ar"
            self.mass = 37211 #MeV
            self.charge = 18
            self.lifetime = None
        elif "potassium" in particleType.lower() or particleType=="K":
            if "40" in particleType:
                self.type = "K-40"
                self.mass = 37226.23 #MeV
                self.charge = 19
                self.lifetime = None
            else:
                self.type = "K-39"
                self.mass = 36294.46 #MeV
                self.charge = 19
                self.lifetime = None
        elif "calcium" in particleType.lower() or particleType=="Ca":
            self.type = "Ca"
            self.mass = 37224.92 #MeV
            self.charge = 20
            self.lifetime = None
        else:
            raise ParticleError("Unrecognized particle type "+str(particleType))


    def buildMomentumFromScalar(self,scalar,kind,otherArgs):
        """Sets particle momentum based on scalar value of momentum, velocity, or energy"""
        if kind=="momentum":
            momentum = scalar
        elif kind=="velocity":
            momentum = scalar*self.mass/c
        elif kind=="energy" or kind=="kinetic":
            if kind=="kinetic":
                scalar += self.mass
            if scalar<self.mass:
                raise ParticleError("Particle energy less than mass ("\
                                    +str(scalar)+"<"+str(self.mass)+")")
            momentum = sqrt(scalar**2 - self.mass**2)
        if ("theta" in otherArgs.keys()) and ("phi" in otherArgs.keys()):
            theta = otherArgs["theta"]
            phi = otherArgs["phi"]
        elif ("zenith" in otherArgs.keys()) and ("azimuth" in otherArgs.keys()):
            theta = pi - otherArgs["zenith"]
            phi = otherArgs["azimuth"] - pi
        else:
            if not("energy" in kind):
                kind = "magnitude of "+kind
            raise ParticleError("Must define angle if "+kind+" is given")
        self.momentum[0] = momentum*sin(theta)*cos(phi)
        self.momentum[1] = momentum*sin(theta)*sin(phi)
        self.momentum[2] = momentum*cos(theta)


    def __getattr__(self,name):
        if name=="Pmag":
            return sqrt(self.momentum[0]**2+self.momentum[1]**2+self.momentum[2]**2)
        if name=="direction" or name=="dir":
            if self.Pmag==0:
                return [0,0,0]
            return [p/self.Pmag for p in self.momentum]
        if name=="theta":
            if self.Pmag==0:
                return 0
            return arccos(self.momentum[2]/self.Pmag)
        if name=="phi":
            if self.momentum[0]==0:
                return 0
            return arctan(self.momentum[1]/self.momentum[0])
        if name=="zenith":
            return pi-self.theta
        if name=="azimuth":
            az = pi+self.phi
            while az>2*pi:
                az -= 2*pi
            return az
        if name=="energy":
            return sqrt(self.mass**2+self.Pmag**2)
        if name=="kinetic" or name=="ke":
            return self.energy-self.mass
        if name=="beta":
            return self.Pmag/self.energy

    # Since __getattr__ is defined, must define these functions to be pickle-able
    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__.update(d)

