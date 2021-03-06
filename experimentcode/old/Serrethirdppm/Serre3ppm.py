# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.11
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Serre3ppm', [dirname(__file__)])
        except ImportError:
            import _Serre3ppm
            return _Serre3ppm
        if fp is not None:
            try:
                _mod = imp.load_module('_Serre3ppm', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Serre3ppm = swig_import_helper()
    del swig_import_helper
else:
    import _Serre3ppm
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def conc(*args):
  return _Serre3ppm.conc(*args)
conc = _Serre3ppm.conc

def mallocPy(*args):
  return _Serre3ppm.mallocPy(*args)
mallocPy = _Serre3ppm.mallocPy

def writetomem(*args):
  return _Serre3ppm.writetomem(*args)
writetomem = _Serre3ppm.writetomem

def readfrommem(*args):
  return _Serre3ppm.readfrommem(*args)
readfrommem = _Serre3ppm.readfrommem

def deallocPy(*args):
  return _Serre3ppm.deallocPy(*args)
deallocPy = _Serre3ppm.deallocPy

def TDMA(*args):
  return _Serre3ppm.TDMA(*args)
TDMA = _Serre3ppm.TDMA

def PENT(*args):
  return _Serre3ppm.PENT(*args)
PENT = _Serre3ppm.PENT

def midpt2ca(*args):
  return _Serre3ppm.midpt2ca(*args)
midpt2ca = _Serre3ppm.midpt2ca

def ca2midpt(*args):
  return _Serre3ppm.ca2midpt(*args)
ca2midpt = _Serre3ppm.ca2midpt

def ufromGh(*args):
  return _Serre3ppm.ufromGh(*args)
ufromGh = _Serre3ppm.ufromGh

def Gfromuh(*args):
  return _Serre3ppm.Gfromuh(*args)
Gfromuh = _Serre3ppm.Gfromuh

def reconstructppm(*args):
  return _Serre3ppm.reconstructppm(*args)
reconstructppm = _Serre3ppm.reconstructppm

def weightsum(*args):
  return _Serre3ppm.weightsum(*args)
weightsum = _Serre3ppm.weightsum

def evolve(*args):
  return _Serre3ppm.evolve(*args)
evolve = _Serre3ppm.evolve

def evolvewrap(*args):
  return _Serre3ppm.evolvewrap(*args)
evolvewrap = _Serre3ppm.evolvewrap
# This file is compatible with both classic and new-style classes.


