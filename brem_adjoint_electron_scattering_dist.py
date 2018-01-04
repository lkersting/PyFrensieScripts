#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.Utility.Interpolation as Interpolation
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
import numpy
import matplotlib.pyplot as plt

Utility.initFrensiePrng()

#datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/frensie/src/packages/test_files/'

source = Teuchos.FileInputSource( datadir + '/cross_sections.xml' )
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList( xml_obj )
data_list = cs_list.get( 'H-Native' )

# -------------------------------------------------------------------------- ##
#  Adjoint Brem Data
# -------------------------------------------------------------------------- ##
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()

###
###  Brem Distribution/Reaction Unit Test Check
###
adjoint_brem_cs = adjoint_data.getAdjointBremsstrahlungElectronCrossSection()

energy = 1e-5
print "energy = ", energy
print '\tcs = ','%.16e' % adjoint_brem_cs[0]

energy = 1e-3
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] <= energy:
        index = i

energy_0 = adjoint_energy_grid[index]
cs_0 = adjoint_brem_cs[index]
energy_1 = adjoint_energy_grid[index+1]
cs_1 = adjoint_brem_cs[index+1]

print "energy = ", energy
cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
print '\tcs = ','%.16e' % cs


energy = 2e-2
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] <= energy:
        index = i

energy_0 = adjoint_energy_grid[index]
cs_0 = adjoint_brem_cs[index]
energy_1 = adjoint_energy_grid[index+1]
cs_1 = adjoint_brem_cs[index+1]

print "energy = ", energy
cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
print '\tcs = ','%.16e' % cs


energy = 20.0
print "energy = ", energy
print '\tcs = ','%.16e' % adjoint_brem_cs[adjoint_brem_cs.size -1]


brem_dist = Collision.createLogLogLogUnitBaseCorrelatedBremsstrahlungDistribution(adjoint_data, 1e-7)

E_in = [1e-6, 1e-5, 1.1e-5, 20.0, 21.0]
E_out = [2.0e-5, 20.2, 1.0, 20.000000201, 22.0]
print "\nEvaluate[ E_in, E_out]"
for i in range(0,len(E_in)):
  pdf = brem_dist.evaluate( E_in[i], E_out[i] )
  print "\teval[",E_in[i],",",E_out[i],"] =\t",'%.16e' % pdf

print "\nEvaluate PDF[ E_in, E_out]"
for i in range(0,len(E_in)):
  pdf = brem_dist.evaluatePDF( E_in[i], E_out[i] )
  print "\teval[",E_in[i],",",E_out[i],"] =\t",'%.16e' % pdf

E_in = [1e-6, 1e-5, 1.1e-5, 20.0, 21.0]
E_out = [2.0e-5, 10.1000050505, 1.0, 20.1000000505, 22.0]
print "\nEvaluate CDF[ E_in, E_out]"
for i in range(0,len(E_in)):
  cdf = brem_dist.evaluateCDF( E_in[i], E_out[i] )
  print "\teval[",E_in[i],",",E_out[i],"] =\t",'%.16e' % cdf

random_numbers = [ 0.5, 0.5 ]
Prng.RandomNumberGenerator.setFakeStream(random_numbers)

incoming_energy = 1.1e-5
outgoing_energy, scattering_angle = brem_dist.sample( incoming_energy )
print "\noutgoing_energy = ",'%.16e' % outgoing_energy
print "scattering_angle = ",'%.16e' % scattering_angle

