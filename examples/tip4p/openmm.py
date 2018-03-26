from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

model = 'tip4p-2005'
resorces = {1:'CPU', 2:'CUDA', 3:'OpenCL'}
choice = 1
temp = 298*unit.kelvin
press = 1.0*unit.atmosphere
dt = 1.0*unit.femtosecond
rc = 9.0*unit.angstrom
nsteps = 10000

pdb = app.PDBFile(model+'.pdb')
forcefield = app.ForceField(model+'.xml')

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=rc,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005)

integrator = mm.LangevinIntegrator(temp, 1.0/unit.picoseconds, dt)
integrator.setConstraintTolerance(0.00001)
system.addForce(mm.MonteCarloBarostat(press, temp, 25))

platform = mm.Platform.getPlatformByName(resorces[choice])
properties = {'Precision': 'mixed'} if choice in (2, 3) else {}

simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

# simulation.reporters.append(app.DCDReporter('trajectory.dcd', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=nsteps, separator='\t'))

simulation.context.setVelocitiesToTemperature(temp)

print('Running Production...')
simulation.step(nsteps)
