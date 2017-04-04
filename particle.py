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
                    momentum = val
                    if ("theta" in kwargs.keys()) and ("phi" in kwargs.keys()):
                        theta = kwargs["theta"]
                        phi = kwargs["phi"]
                    elif ("zenith" in kwargs.keys()) and ("azimuth" in kwargs.keys()):
                        theta = pi - kwargs["zenith"]
                        phi = kwargs["azimuth"] - pi
                    else:
                        raise ParticleError("Must define angle if magnitude of momentum is given")
                    self.momentum[0] = momentum*sin(theta)*cos(phi)
                    self.momentum[1] = momentum*sin(theta)*sin(phi)
                    self.momentum[2] = momentum*cos(theta)
            if "vel" in key:
                if isinstance(val,list) and len(val)==3:
                    self.momentum = [v*self.mass/c for v in val]
                elif isinstance(val,float) or isinstance(val,int):
                    momentum = val*self.mass/c
                    if ("theta" in kwargs.keys()) and ("phi" in kwargs.keys()):
                        theta = kwargs["theta"]
                        phi = kwargs["phi"]
                    elif ("zenith" in kwargs.keys()) and ("azimuth" in kwargs.keys()):
                        theta = pi - kwargs["zenith"]
                        phi = kwargs["azimuth"] - pi
                    else:
                        raise ParticleError("Must define angle if magnitude of velocity is given")
                    self.momentum[0] = momentum*sin(theta)*cos(phi)
                    self.momentum[1] = momentum*sin(theta)*sin(phi)
                    self.momentum[2] = momentum*cos(theta)

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


    def __getattr__(self,name):
        if name=="theta":
            p_mag = sqrt(self.momentum[0]**2+self.momentum[1]**2+self.momentum[2]**2)
            if p_mag==0:
                return 0
            return arccos(self.momentum[2]/p_mag)
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

