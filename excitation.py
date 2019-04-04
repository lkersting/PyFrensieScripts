#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.MonteCarlo.Collision as Collision
import PyFrensie.MonteCarlo.Electron as Electron
import numpy
import matplotlib.pyplot as plt

Utility.initFrensiePrng()

### -------------------------------------------------------------------------- ##
###  Forward Atomic Excitation
### -------------------------------------------------------------------------- ##
native_file_name = '/home/software/mcnpdata/native/epr/epr_native_1_v1.xml'
# native_file_name = '/home/lkersting/frensie/src/packages/test_files/native/test_epr_1_native.xml'
native_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )

dist = Electron.createAtomicExcitationDistribution( native_data )

energies = [0.01, 1.0, 5.0, 9.50367370e-03]

# print "\t -- Evaluate --"
# for energy in energies:
#   print "Energy = ",energy
#   if energy == 1e5:
#       outgoing_energies = [1e-5, 1.0, 10.0]
#   else:
#       outgoing_energies = [1e-5, 1e-4, 5e-4]
#   for e_out in outgoing_energies:
#     pdf = dist.evaluate( energy, e_out )
#     print '\teval[','%.6e' %e_out,']\t= ','%.16e' % pdf

# print "\n\t-- Evaluate PDF --"
# for energy in energies:
#   print "Energy = ",energy
#   if energy == 1e5:
#       outgoing_energies = [1e-5, 1.0, 10.0]
#   else:
#       outgoing_energies = [1e-5, 1e-4, 5e-4]
#   for e_out in outgoing_energies:
#     pdf = dist.evaluatePDF( energy, e_out )
#     print '\t pdf[','%.6e' %e_out,']\t= ','%.16e' % pdf


# print "\n\t-- Evaluate CDF --"
# for energy in energies:
#   print "Energy = ",energy
#   if energy == 1e5:
#       outgoing_energies = [1e-5, 1.0, 10.0]
#   else:
#       outgoing_energies = [1e-5, 1e-4, 5e-4]
#   for e_out in outgoing_energies:
#     scattering_e_out_cosine = -0.01
#     cdf = dist.evaluateCDF( energy, e_out )
#     print '\t cdf[','%.6e' %e_out,']\t= ','%.16e' % cdf

N = 20
random_numbers = [ 0.5 ]*N
Prng.RandomNumberGenerator.setFakeStream(random_numbers)

for energy in energies:

  print "\n\t-- Sample at energy", energy, "--"
  e_in = energy
  for i in range(len(random_numbers)):
      e_out,angle = dist.sample( e_in )
      print '%.20e' % e_out
      e_in = e_out