"""Code for particle class"""
from constants import c, pi
from numpy import sqrt,sin,cos,arccos,arctan

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
            if "mom" in key:
                if isinstance(val,list) and len(val)==3:
                    self.momentum = [float(p) for p in val]
                elif isinstance(val,float) or isinstance(val,int):
                    self.buildMomentumFromScalar(val,"momentum",kwargs)
            if "vel" in key:
                if isinstance(val,list) and len(val)==3:
                    self.momentum = [v*self.mass/c for v in val]
                elif isinstance(val,float) or isinstance(val,int):
                    self.buildMomentumFromScalar(val,"velocity",kwargs)
            if key=="energy" or key=="E":
                self.buildMomentumFromScalar(val,"energy",kwargs)


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
        elif particleType.lower()=="nitrogen" or particleType=="N":
            self.type = "N"
            self.mass = 13047.2 #MeV
            self.charge = 7
        elif particleType.lower()=="oxygen" or particleType=="O":
            self.type = "O"
            self.mass = 14903.3 #MeV
            self.charge = 8
        elif particleType.lower()=="argon" or particleType=="Ar":
            self.type = "Ar"
            self.mass = 37211 #MeV
            self.charge = 18
        else:
            raise ParticleError("Unrecognized particle type "+str(particleType))


    def buildMomentumFromScalar(self,scalar,kind,otherArgs):
        """Sets particle momentum based on scalar value of momentum, velocity, or energy"""
        if kind=="momentum":
            momentum = scalar
        elif kind=="velocity":
            momentum = scalar*self.mass/c
        elif kind=="energy":
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

