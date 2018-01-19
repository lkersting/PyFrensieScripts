#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
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

# -------------------------------------------------------------------------- ##
#  Brem Data
# -------------------------------------------------------------------------- ##

# Possible elements ['Al-Native', 'H-Native', 'Pb-Native']
elements = ['Al-Native']
# Possible Interpolation Schemes ["Unit-base", "Unit-base Correlated", "Correlated"]
schemes = ["Unit-base Correlated"]
# Possible interpolation schemes ["LogLogLog", "LogLogLog", "LinLinLin"]
interps = ["LogLogLog"]
# Possible energies [1e-5, 1e-3, 1e5 ]
energies = [1e-5, 1e-3, 1e5 ]
# Only do the following shell
some_shells_only = True
selected_shells = [6]
# Evaluation tolerance
tol = 1e-12

for z in elements:
  print "\n----------------------------"
  print "-----", z, "Tests -----"
  print "----------------------------"
  data_list = cs_list.get(z)
  file_name = datadir + data_list.get('electroatomic_file_path')
  print file_name
  native_data = Native.ElectronPhotonRelaxationDataContainer(file_name)

  print "----------------------------\n"

  for interp in interps:
    for scheme in schemes:
      print "\n-----", interp, "-", scheme, "-----\n"
      if interp == "LogLogLog":
        if scheme == "Unit-base":
          dist = Collision.createLogLogLogUnitBaseBremsstrahlungDistribution(native_data, tol)
        if scheme == "Unit-base Correlated":
          dist = Collision.createLogLogLogUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
        if scheme == "Correlated":
          dist = Collision.createLogLogLogCorrelatedBremsstrahlungDistribution(native_data, tol)
      elif interp == "LinLinLin":
        if scheme == "Unit-base":
          dist = Collision.createLinLinLinUnitBaseBremsstrahlungDistribution(native_data, tol)
        elif scheme == "Unit-base Correlated":
          dist = Collision.createLinLinLinUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
        elif scheme == "Correlated":
          dist = Collision.createLinLinLinCorrelatedBremsstrahlungDistribution(native_data, tol)
      elif interp == "LinLinLog":
        if scheme == "Unit-base":
          dist = Collision.createLinLinLogUnitBaseBremsstrahlungDistribution(native_data, tol)
        elif scheme == "Unit-base Correlated":
          dist = Collision.createLinLinLogUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
        elif scheme == "Correlated":
          dist = Collision.createLinLinLogCorrelatedBremsstrahlungDistribution(native_data, tol)

      incoming_energies = [0.02, 9.0e-4, 1.0e5]
      outgoing_energies = [1e-7, 9.0e-4, 2e4]
      print "\t -- Evaluate --"
      print "eval[e_in, e_out]"
      for i in range(0, len(incoming_energies)):
        e_in = incoming_energies[i]
        e_out = outgoing_energies[i]
        pdf = dist.evaluate(e_in, e_out)
        print '\t', i+1, ') eval[', '%.6e' %e_in, '%.6e' %e_out, ']\t= ', '%.16e' % pdf

      print "\n\t-- Evaluate PDF --"
      print "pdf[e_in, e_out]"
      for i in range(0, len(incoming_energies)):
        e_in = incoming_energies[i]
        e_out = outgoing_energies[i]
        pdf = dist.evaluatePDF(e_in, e_out)
        print '\t', i+1, ') pdf[', '%.6e' %e_in, '%.6e' %e_out, ']\t= ', '%.16e' % pdf

      print "\n\t-- Evaluate CDF --"
      print "cdf[e_in, e_out]"
      for i in range(0, len(incoming_energies)):
        e_in = incoming_energies[i]
        e_out = outgoing_energies[i]
        cdf = dist.evaluateCDF(e_in, e_out)
        print '\t', i+1, ') cdf[', '%.6e' %e_in, '%.6e' %e_out, ']\t= ', '%.16e' % cdf


      energies = [6.041e-05, 1e-3, 1e-4]
      for energy in energies:
        if energy == 1e-3:
            random_numbers = [0.0, 0.0, 1.0 - 1e-15, 1.0 - 1e-15]
        elif energy == 6.041e-05:
          random_numbers = [0.0, 0.0, 1.0 - 1e-15, 0.0]
        else:
            random_numbers = [0.5, 0.5]

        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

        print "\n\t-- Sample --"
        print "Energy = ", energy
        for i in range(0, len(random_numbers)/2):
            e_out, angle = dist.sample(energy)
            print '\te_out = ', '%.16e' % e_out, '\tangle = ', '%.16e' % angle

energy_grid = native_data.getBremsstrahlungEnergyGrid()
linlinlin_brem_dist = Collision.createLinLinLinUnitBaseCorrelatedBremsstrahlungDistribution(native_data, 1e-7)
linlinlog_brem_dist = Collision.createLinLinLogUnitBaseCorrelatedBremsstrahlungDistribution(native_data, 1e-7)
logloglog_brem_dist = Collision.createLogLogLogUnitBaseCorrelatedBremsstrahlungDistribution(native_data, 1e-7)

energy = 0.0009

index = 0
for i in range(0, energy_grid.size ):
    if energy_grid[i] <= energy:
        index = i

energy_0 = energy_grid[index]
energy_1 = energy_grid[index+1]
print energy_0
print native_data.getBremsstrahlungPhotonEnergy( energy_0 )
print energy_1
print native_data.getBremsstrahlungPhotonEnergy( energy_1 )

cdf_value = [0.5]
for cdf in cdf_value:
    random_numbers = [ cdf, cdf, cdf, cdf, cdf, cdf, cdf ]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    print "\n\t--- Lower Incoming Electron Energy ",energy_0," ---"
    outgoing_E_0, angle_0 = linlinlin_brem_dist.sample( energy_0 )
    linlinlog_outgoing_E_0, linlinlog_angle_0 = linlinlog_brem_dist.sample( energy_0 )
    print "\tLin-Lin-Lin: E_0 = ",'%.18e' % outgoing_E_0,"\tangle_0 = ",'%.18e' % angle_0
    print "\tLin-Lin-Log: E_0 = ",'%.18e' % linlinlog_outgoing_E_0,"\tangle_0 = ",'%.18e' % linlinlog_angle_0

    print "\n\t--- Upper Incoming Electron Energy ",energy_1," ---"
    outgoing_E_1, angle_1 = linlinlin_brem_dist.sample( energy_1 )
    linlinlog_outgoing_E_1, linlinlog_angle_1 = linlinlog_brem_dist.sample( energy_1 )
    print "\tLin-Lin-Lin: E_1 = ",'%.18e' % outgoing_E_1,"\tangle_1 = ",'%.18e' % angle_1
    print "\tLin-Lin-Log: E_1 = ",'%.18e' % linlinlog_outgoing_E_1,"\tangle_1 = ",'%.18e' % linlinlog_angle_1


    E_in_log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
    E_in_lin_interp = (energy-energy_0)/(energy_1-energy_0)

    loglog_outgoing_E = numpy.exp(numpy.log(outgoing_E_0) + numpy.log(outgoing_E_1/outgoing_E_0)*E_in_log_interp)
    linlog_outgoing_E = outgoing_E_0 + (outgoing_E_1 - outgoing_E_0)*E_in_log_interp
    linlin_outgoing_E = outgoing_E_0 + (outgoing_E_1 - outgoing_E_0)*E_in_lin_interp

    print "\n\t--- Incoming Electron Energy ",energy," ---"
    print "\tLinLin Interp Outgoing Energy: ",'%.18e' % linlin_outgoing_E
    print "\tLinLog Interp Outgoing Energy: ",'%.18e' % linlog_outgoing_E
    print "\tLogLog Interp Outgoing Energy: ",'%.18e' % loglog_outgoing_E

    sample_linlinlin_E, sample_linlinlin_angle = linlinlin_brem_dist.sample( energy )
    print "\tLin-Lin-Lin: Outgoing Energy: = ",'%.18e' % sample_linlinlin_E,"\tangle = ",'%.18e' % sample_linlinlin_angle

    sample_linlinlog_E, sample_linlinlog_angle = linlinlog_brem_dist.sample( energy )
    print "\tLin-Lin-Log: Outgoing Energy: = ",'%.18e' % sample_linlinlog_E,"\tangle = ",'%.18e' % sample_linlinlog_angle

    sample_logloglog_E, sample_logloglog_angle = logloglog_brem_dist.sample( energy )
    print "\tLog-Log-Log: Outgoing Energy: = ",'%.18e' % sample_logloglog_E,"\tangle = ",'%.18e' % sample_logloglog_angle

