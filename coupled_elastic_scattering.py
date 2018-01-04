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
###  Coupled Distribution/Reaction Unit Test Check
###
interps = ["LinLinLog", "LogLogLog", "LinLinLin"]
methods = ["Simplified Union","One D Union","Two D Union"]
energies = [1.0e+5,6.625E+01,200.0]

interps = ["LogLogLog"]
methods = ["Simplified Union"]
energies = [1.0e+5,6.625E+01,200.0]
for interp in interps:
  print "\n-----", interp ,"-----\n"
  for method in methods:
      coupled_dist = Collision.createLinLinLogCoupledElasticDistribution( native_data, method, 1e-15 )
      if interp == "LogLogLog":
        coupled_dist = Collision.createLogLogLogCorrelatedCoupledElasticDistribution( native_data, method, 1e-15 )
      elif interp == "LinLinLin":
        coupled_dist = Collision.createLinLinLinCorrelatedCoupledElasticDistribution( native_data, method, 1e-15 )


      print "\n----- ",method," -----\n"

      print "\t -- Evaluate --"
      angles = [-0.01, 0.0, 0.71, 0.999999, 1.0]
      for energy in energies:
        print "Energy = ",energy
        for angle in angles:
          pdf = coupled_dist.evaluate( energy, angle )
          print '\teval[','%.6e' %angle,']\t= ','%.16e' % pdf

      print "\n\t-- Evaluate PDF --"
      angles = [-0.01, 0.0, 0.71, 0.999999, 1.0]
      for energy in energies:
        print "Energy = ",energy
        for angle in angles:
          pdf = coupled_dist.evaluatePDF( energy, angle )
          print '\t pdf[','%.6e' %angle,']\t= ','%.16e' % pdf


      print "\n\t-- Evaluate CDF --"
      angles = [-0.01, 0.0, 0.71, 0.999995, 0.999999, 1.0]
      for energy in energies:
        print "Energy = ",energy
        for angle in angles:
          scattering_angle_cosine = -0.01
          cdf = coupled_dist.evaluateCDF( energy, angle )
          print '\t cdf[','%.6e' %angle,']\t= ','%.16e' % cdf


      energies = [66.25, 1e-4, 2e2]
      for energy in energies:
        cutoff_cdf = coupled_dist.evaluateCDF( energy, 0.999999 )
        random_numbers = [ 0.0, 1.0e-3, cutoff_cdf, 2.26e-3, 0.5, 1.0 - 1e-15]
        if energy == 1e-4:
          random_numbers = [ 0.0, 1.0e-3, 0.5, 0.9, 0.999999, 1.0-1e-15]
        elif energy == 2e2:
          random_numbers = [ 0.0, 1e-5, 1.0e-4, 2.510603e-04, cutoff_cdf, 1.0 - 1e-15]
        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

        print "\n\t-- Sample --"
        print "Energy = ",energy, "\tcutoff_cdf = ", '%.16e' % cutoff_cdf
        for random_number in random_numbers:
          energy,angle = coupled_dist.sample( energy )
          print '\tangle[','%.6e' %random_number,'] = ','%.16e' % angle