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
native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
h_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )

data_list = cs_list.get( 'Pb-Native' )
native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
pb_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )

###
###  Construct Coupled Distribution
###
interps = ["LinLinLog", "LogLogLog", "LinLinLin"]
methods = ["Simplified Union","One D Union","Two D Union"]


interps = ["LogLogLog", "LinLinLin"]
methods = ["Two D Union"]
for interp in interps:
  print "\n-----", interp ,"-----"
  for method in methods:
      coupled_dist = Collision.createLinLinLogCorrelatedCoupledElasticDistribution( pb_data, method, 1e-15 )
      if interp == "LogLogLog":
        coupled_dist = Collision.createLogLogLogCorrelatedCoupledElasticDistribution( pb_data, method, 1e-15 )
      elif interp == "LinLinLin":
        coupled_dist = Collision.createLinLinLinCorrelatedCoupledElasticDistribution( pb_data, method, 1e-15 )


      print "\n----- ",method," -----\n"

      print "\t -- evaluateLegendreExpandedRutherford --\n"
      etas = [2.51317958942017e3, 2.68213671998009, 4.14887699806239e-14, 2.51317958942017e3]
      angles = [[0.999999], [1.0, 1.0, 1.0], [0.999999], [1.0]]
      energies = [[1.0e-5], [1.0e-4, 5.5e-4, 1.0e-3], [1.0e5], [1.0e-5]]
      for i in range(0,len(etas)):
        print "Eta = ",etas[i]
        for j in range(0, len(angles[i])):
          cutoff_pdf = coupled_dist.evaluatePDFAtCutoff( energies[i][j] )
          pdf = coupled_dist.evaluateScreenedRutherfordPDF( angles[i][j], etas[i], cutoff_pdf )
          print '\tpdf[','%.6e' %angles[i][j],',','%.6e' %energies[i][j],']\t= ','%.20e' % pdf

      print "\n\t -- evaluateLegendreExpandedPDF --\n"
      energies = [1.0e-5, 1.0e-4, 5.5e-4, 1.0e-3, 1.0e5]
      angles = [[-1.0, 0.999999], [-1.0, 0.999999], [0.999999], [-1.0, 0.999999], [0.9999979]]
      for i in range(0,len(energies)):
        print "Energy = ",energies[i]
        for j in range(0, len(angles[i])):
          pdf = coupled_dist.evaluatePDF( energies[i],angles[i][j] )
          print '\tpdf[','%.6e' %angles[i][j],']\t= ','%.20e' % pdf


      #   print "Energy = ",energy
      #   for angle in angles:
      #     pdf = coupled_dist.evaluateScreenedRutherfordPDF( energy, angle )
      #     print '\teval[','%.6e' %angle,']\t= ','%.16e' % pdf

      # print "\n\t-- Evaluate PDF --"
      # angles = [-0.01, 0.0, 0.71, 0.999999, 1.0]
      # for energy in energies:
      #   print "Energy = ",energy
      #   for angle in angles:
      #     pdf = coupled_dist.evaluateScreenedRutherfordPDF( energy, angle )
      #     print '\t pdf[','%.6e' %angle,']\t= ','%.16e' % pdf


      # print "\n\t-- Evaluate CDF --"
      # angles = [-0.01, 0.0, 0.71, 0.999995, 0.999999, 1.0]
      # for energy in energies:
      #   print "Energy = ",energy
      #   for angle in angles:
      #     scattering_angle_cosine = -0.01
      #     cdf = coupled_dist.evaluateCDF( energy, angle )
      #     print '\t cdf[','%.6e' %angle,']\t= ','%.16e' % cdf


      # energies = [66.25, 1e-4, 2e2]
      # for energy in energies:
      #   cutoff_cdf = coupled_dist.evaluateCDF( energy, 0.999999 )
      #   random_numbers = [ 0.0, 1.0e-3, cutoff_cdf, 2.26e-3, 0.5, 1.0 - 1e-15]
      #   if energy == 1e-4:
      #     random_numbers = [ 0.0, 1.0e-3, 0.5, 0.9, 0.999999, 1.0-1e-15]
      #   elif energy == 2e2:
      #     random_numbers = [ 0.0, 1e-5, 1.0e-4, 2.510603e-04, cutoff_cdf, 1.0 - 1e-15]
      #   Prng.RandomNumberGenerator.setFakeStream(random_numbers)

      #   print "\n\t-- Sample --"
      #   print "Energy = ",energy, "\tcutoff_cdf = ", '%.16e' % cutoff_cdf
      #   for random_number in random_numbers:
      #     energy,angle = coupled_dist.sample( energy )
      #     print '\tangle[','%.6e' %random_number,'] = ','%.16e' % angle