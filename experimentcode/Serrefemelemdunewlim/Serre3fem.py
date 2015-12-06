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
            fp, pathname, description = imp.find_module('_Serre3fem', [dirname(__file__)])
        except ImportError:
            import _Serre3fem
            return _Serre3fem
        if fp is not None:
            try:
                _mod = imp.load_module('_Serre3fem', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Serre3fem = swig_import_helper()
    del swig_import_helper
else:
    import _Serre3fem
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
  return _Serre3fem.conc(*args)
conc = _Serre3fem.conc

def mallocPy(*args):
  return _Serre3fem.mallocPy(*args)
mallocPy = _Serre3fem.mallocPy

def writetomem(*args):
  return _Serre3fem.writetomem(*args)
writetomem = _Serre3fem.writetomem

def readfrommem(*args):
  return _Serre3fem.readfrommem(*args)
readfrommem = _Serre3fem.readfrommem

def deallocPy(*args):
  return _Serre3fem.deallocPy(*args)
deallocPy = _Serre3fem.deallocPy

def TDMA(*args):
  return _Serre3fem.TDMA(*args)
TDMA = _Serre3fem.TDMA

def PENT(*args):
  return _Serre3fem.PENT(*args)
PENT = _Serre3fem.PENT

def midpt2ca(*args):
  return _Serre3fem.midpt2ca(*args)
midpt2ca = _Serre3fem.midpt2ca

def ca2midpt(*args):
  return _Serre3fem.ca2midpt(*args)
ca2midpt = _Serre3fem.ca2midpt

def ufromGh(*args):
  return _Serre3fem.ufromGh(*args)
ufromGh = _Serre3fem.ufromGh

def Gfromuh(*args):
  return _Serre3fem.Gfromuh(*args)
Gfromuh = _Serre3fem.Gfromuh

def weightsum(*args):
  return _Serre3fem.weightsum(*args)
weightsum = _Serre3fem.weightsum

def evolve(*args):
  return _Serre3fem.evolve(*args)
evolve = _Serre3fem.evolve

def evolvewrap(*args):
  return _Serre3fem.evolvewrap(*args)
evolvewrap = _Serre3fem.evolvewrap
# This file is compatible with both classic and new-style classes.


