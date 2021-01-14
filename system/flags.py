#/usr/bin/env python


### Flags for systems in SIM suite.
### coded by MSS



#algorithm-maintained flags
FlagDflts = {"CalcForce" : False,
             "CalcVirial" : False,
             "CalcDUParam" : False,
             "CalcDWParam" : False,
             "ConservesMomentum" : False,
             "CalcFluct" : False,}
                
#flags to turn off for minimal calculations
OffFlags = ["CalcDUParam", "CalcDWParam", "CalcForce", "CalcVirial", "CalcFluct"]



class FlagsClass(object):

    def __init__(self):
        for (k,v) in FlagDflts.iteritems():
            setattr(self, k, v)
            
    def CalcsOff(self):
        for f in OffFlags:
            if hasattr(self, f):
                setattr(self, f, False)
            else:
                raise ValueError("Don't recognize flag %s" % f)
                
    
        
            


