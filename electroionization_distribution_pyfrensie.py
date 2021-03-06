#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.Utility.Interpolation as Interpolation
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
import numpy

Utility.initFrensiePrng()

#datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/frensie/src/packages/test_files/'

source = Teuchos.FileInputSource( datadir + '/cross_sections.xml' )
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList( xml_obj )

# -------------------------------------------------------------------------- ##
#  Electroionization Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'H-Native' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
energy_grid = native_data.getElectronEnergyGrid()

###
###  Electroionization Reaction Unit Test Check
###
print "\n----- H -----"
subshells = native_data.getSubshells()
shell = subshells[0]
binding_energy = native_data.getSubshellBindingEnergy( shell )

interps = ["LinLinLog", "LogLogLog", "LinLinLin"]
energies = [1.0e+5, 1e-3]

interps = ["LogLogLog"]
energies = [1.0e+5, 1e-3]
for interp in interps:
    print "\n-----", interp ,"-----\n"
    ionization_dist = Collision.createLinLinLogUnitBaseCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-15)
    if interp == "LogLogLog":
      ionization_dist = Collision.createLogLogLogUnitBaseCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-15)
    elif interp == "LinLinLin":
      ionization_dist = Collision.createLinLinLinUnitBaseCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-15)


    print "\nshell = ", shell

    print "\t -- Evaluate --"
    for energy in energies:
      print "Energy = ",energy
      if energy == 1e5:
          outgoing_energies = [1e-5, 1.0, 10.0]
      else:
          outgoing_energies = [1e-5, 1e-4, 5e-4]
      for e_out in outgoing_energies:
        pdf = ionization_dist.evaluate( energy, e_out )
        print '\teval[','%.6e' %e_out,']\t= ','%.16e' % pdf


    print "\n\t-- Evaluate PDF --"
    for energy in energies:
      print "Energy = ",energy
      if energy == 1e5:
          outgoing_energies = [1e-5, 1.0, 10.0]
      else:
          outgoing_energies = [1e-5, 1e-4, 5e-4]
      for e_out in outgoing_energies:
        pdf = ionization_dist.evaluatePDF( energy, e_out )
        print '\t pdf[','%.6e' %e_out,']\t= ','%.16e' % pdf


    print "\n\t-- Evaluate CDF --"
    for energy in energies:
      print "Energy = ",energy
      if energy == 1e5:
          outgoing_energies = [1e-5, 1.0, 10.0]
      else:
          outgoing_energies = [1e-5, 1e-4, 5e-4]
      for e_out in outgoing_energies:
        scattering_e_out_cosine = -0.01
        cdf = ionization_dist.evaluateCDF( energy, e_out )
        print '\t cdf[','%.6e' %e_out,']\t= ','%.16e' % cdf


    energies = [1e-3, 1e5]
    for energy in energies:
      random_numbers = [ 0.0, 1.0 - 1e-15]

      Prng.RandomNumberGenerator.setFakeStream(random_numbers)

      print "\n\t-- Sample --"
      print "Energy = ",energy
      for i in range(0, len(random_numbers)):
          e_out,angle = ionization_dist.sample( energy )
          print '\te_out = ','%.16e' % e_out,'\tangle = ','%.16e' % angle