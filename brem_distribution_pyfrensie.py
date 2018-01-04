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

### -------------------------------------------------------------------------- ##
###  Forward Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )
energy_grid = native_data.getElectronEnergyGrid()

tot_elastic_cs = native_data.getTotalElasticCrossSection()
cutoff_cs = native_data.getCutoffElasticCrossSection()
screen_rutherford_cs = native_data.getScreenedRutherfordElasticCrossSection()
screen_rutherford_index = native_data.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()


###
###  Brem Distribution/Reaction Unit Test Check
###
interps = ["LinLinLog", "LogLogLog", "LinLinLin"]
energies = [1.0e+5, 1e-3]

interps = ["LogLogLog"]
energies = [1.0e+5, 1e-3]
for interp in interps:
    print "\n-----", interp ,"-----\n"
    brem_dist = Collision.createLinLinLogUnitBaseCorrelatedBremsstrahlungDistribution( native_data, 1e-15 )
    if interp == "LogLogLog":
      brem_dist = Collision.createLogLogLogUnitBaseCorrelatedBremsstrahlungDistribution( native_data, 1, 1e-15 )
    elif interp == "LinLinLin":
      brem_dist = Collision.createLinLinLinUnitBaseCorrelatedBremsstrahlungDistribution( native_data, 1e-15 )

    print "\t -- Evaluate --"
    for energy in energies:
      print "Energy = ",energy
      if energy == 1e5:
          outgoing_energies = [1e-5, 1.0, 10.0]
      else:
          outgoing_energies = [1e-5, 1e-4, 5e-4]
      for e_out in outgoing_energies:
        pdf = brem_dist.evaluate( energy, e_out )
        print '\teval[','%.6e' %e_out,']\t= ','%.16e' % pdf

    print "\n\t-- Evaluate PDF --"
    for energy in energies:
      print "Energy = ",energy
      if energy == 1e5:
          outgoing_energies = [1e-5, 1.0, 10.0]
      else:
          outgoing_energies = [1e-5, 1e-4, 5e-4]
      for e_out in outgoing_energies:
        pdf = brem_dist.evaluatePDF( energy, e_out )
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
        cdf = brem_dist.evaluateCDF( energy, e_out )
        print '\t cdf[','%.6e' %e_out,']\t= ','%.16e' % cdf


    energies = [1e-3, 1e-4]
    for energy in energies:
      if energy == 1e-3:
          random_numbers = [ 0.0, 0.0, 1.0 - 1e-15, 1.0 - 1e-15]
      else:
          random_numbers = [ 0.5, 0.5]

      Prng.RandomNumberGenerator.setFakeStream(random_numbers)

      print "\n\t-- Sample --"
      print "Energy = ",energy
      for i in range(0, len(random_numbers)/2):
          e_out,angle = brem_dist.sample( energy )
          print '\te_out = ','%.16e' % e_out,'\tangle = ','%.16e' % angle