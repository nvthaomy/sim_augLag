#/usr/bin/env python


### Standard integration routine objects in SIM suite.
### coded by MSS


import velverlet
import montecarlo
import mcmoves
import mcadddel
import mcscalevol
import mcsemigrand


def Add(Int):
    """Adds standard integrators to system object Sys."""
    Sys = Int.Sys
    TimeStep = 0.001 * Sys.Units.TimeScale
    Method = Int.AddMethod(velverlet.VVQuench, TimeStep = TimeStep)
    Method = Int.AddMethod(velverlet.VVIntegrate, TimeStep = TimeStep)
    Method.Thermostat = Method.ThermostatAndersen
    Method.AndersenStepFreq = 100
    Moves = [mcmoves.CanonicalDisplacements(Sys), 
             mcadddel.AddDeleteMolecules(Sys),
             mcscalevol.ScaleVolume(Sys),
             mcsemigrand.MutateMolecules(Sys)]
    Method = Int.AddMethod(montecarlo.MonteCarlo, Moves = Moves)  
    #set the default
    Int.Method = Int.Methods.VVIntegrate