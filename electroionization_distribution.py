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

# Possible elements ['Al-Native', 'H-Native', 'Pb-Native']
elements = ['Al-Native']
# Possible Interpolation Schemes ["Unit-base", "Unit-base Correlated", "Correlated"]
schemes = ["Unit-base Correlated", "Unit-base"]
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
  data_list = cs_list.get( z )
  file_name = datadir + data_list.get( 'electroatomic_file_path' )
  print file_name
  native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
  subshells = native_data.getSubshells()

  if some_shells_only:
    shells = selected_shells
  else:
    shells = subshells

  for shell in shells:

    binding_energy = native_data.getSubshellBindingEnergy( shell )
    print "\nshell = ", shell, "\tbinding energy =", binding_energy
    print "----------------------------\n"

    for interp in interps:
      for scheme in schemes:
        print "\n-----", interp, "-", scheme, "-----\n"
        if interp == "LogLogLog":
          if scheme == "Unit-base":
            dist = Collision.createLogLogLogUnitBaseElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
          if scheme == "Unit-base Correlated":
            dist = Collision.createLogLogLogUnitBaseCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
          if scheme == "Correlated":
            dist = Collision.createLogLogLogCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
        elif interp == "LinLinLin":
          if scheme == "Unit-base":
            dist = Collision.createLinLinLinUnitBaseElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
          elif scheme == "Unit-base Correlated":
            dist = Collision.createLinLinLinUnitBaseCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
          elif scheme == "Correlated":
            dist = Collision.createLinLinLinCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
        elif interp == "LinLinLog":
          if scheme == "Unit-base":
            dist = Collision.createLinLinLogUnitBaseElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
          elif scheme == "Unit-base Correlated":
            dist = Collision.createLinLinLogUnitBaseCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
          elif scheme == "Correlated":
            dist = Collision.createLinLinLogCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)

        incoming_energies = [8.829e-2, 9.12175e-2, 1e-1, 1.0, 1.0, 1e5, 1e5]
        outgoing_energies = [1e-7, 4.275e-4, 1e-2, 1.33136131511529e-1, 9.7163E-02, 1.752970e2, 5e4]
        print "\t -- Evaluate --"
        print "eval[e_in, e_out]"
        for i in range(0, len(incoming_energies)):
          e_in = incoming_energies[i]
          e_out = outgoing_energies[i]
          pdf = dist.evaluate(e_in, e_out)
          print '\t', i, 'eval[', '%.6e' %e_in, '%.6e' %e_out, ']\t= ', '%.16e' % pdf

        print "\n\t-- Evaluate PDF --"
        print "pdf[e_in, e_out]"
        for i in range(0, len(incoming_energies)):
          e_in = incoming_energies[i]
          e_out = outgoing_energies[i]
          pdf = dist.evaluatePDF(e_in, e_out)
          print '\t', i, 'pdf[', '%.6e' %e_in, '%.6e' %e_out, ']\t= ', '%.16e' % pdf

        print "\n\t-- Evaluate CDF --"
        print "cdf[e_in, e_out]"
        for i in range(0, len(incoming_energies)):
          e_in = incoming_energies[i]
          e_out = outgoing_energies[i]
          cdf = dist.evaluateCDF(e_in, e_out)
          print '\t', i, 'cdf[', '%.6e' %e_in, '%.6e' %e_out, ']\t= ', '%.16e' % cdf


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


    #     ionization_dist = Collision.createLogLogLogUnitBaseElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-15)

    #     print "\n--- sample ---";

    #     energies = [1.5, 15.7, binding_energy + 1e-10, 5e4]

    #     e_out = 0.0
    #     knock_on_energy = 0.0
    #     scattering_angle_cosine = 0.0
    #     knock_on_angle_cosine = 0.0

    #     for i in range(0,len(energies)):

    #        random_numbers = [0.0, 1.0-1e-12]
    #        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    #        print "\n\tenergies[i] =",energies[i],"\trandom_number =",random_numbers[0]
    #        knock_on_energy, knock_on_angle_cosine = ionization_dist.sample( energies[i] )
    #        print "knock-on mu   =",'%.16e' % knock_on_angle_cosine, "\tknock-on energy =",'%.16e' % knock_on_energy


    #        print "\tenergies[i] =",energies[i],"\trandom_number =",random_numbers[1]
    #        knock_on_energy, knock_on_angle_cosine = ionization_dist.sample( energies[i] )
    #        print "knock-on mu   =",'%.16e' % knock_on_angle_cosine, "\tknock-on energy =",'%.16e' % knock_on_energy


    #        random_numbers = [0.0, 1.0-1e-12]
    #        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    #        e_out, knock_on_energy, scattering_angle_cosine, knock_on_angle_cosine = ionization_dist.samplePrimaryAndSecondary( energies[i] )
    #        print "\tenergies[i] =",energies[i],"\trandom_number =",random_numbers[0]
    #        print "knock-on mu   =",'%.16e' % knock_on_angle_cosine, "\tknock-on energy =",'%.16e' % knock_on_energy
    #        print "scattering mu =",'%.16e' % scattering_angle_cosine, "\toutgoing energy =",'%.16e' % e_out


    #        e_out, knock_on_energy, scattering_angle_cosine, knock_on_angle_cosine = ionization_dist.samplePrimaryAndSecondary( energies[i] )
    #        print "\tenergies[i] =",energies[i],"\trandom_number =",random_numbers[1]
    #        print "knock-on mu   =",'%.16e' % knock_on_angle_cosine, "\tknock-on energy =",'%.16e' % knock_on_energy
    #        print "scattering mu =",'%.16e' % scattering_angle_cosine, "\toutgoing energy =",'%.16e' % e_out




    #     energy_grid = native_data.getElectroionizationEnergyGrid(shell)
    #     linlinlin_dist = Collision.createLinLinLinCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-15)
    #     linlinlog_dist = Collision.createLinLinLogUnitBaseCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-15)
    #     logloglog_dist = Collision.createLogLogLogCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-15)

    #     energy = 6.041e-05

    #     index = 0
    #     for i in range(0, energy_grid.size ):
    #         if energy_grid[i] <= energy:
    #             index = i

    #     energy_0 = energy_grid[index]
    #     energy_1 = energy_grid[index+1]
    #     print energy_0
    #     recoil_0 = native_data.getElectroionizationRecoilEnergy( shell, energy_0 )
    #     print recoil_0
    #     print energy_1
    #     recoil_1 = native_data.getElectroionizationRecoilEnergy( shell, energy_1 )
    #     print recoil_1

    #     cdf_value = [1.0-1e-15]
    #     for cdf in cdf_value:
    #         random_numbers = [ cdf, cdf, cdf, cdf, cdf, cdf, cdf ]
    #         Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    #         print "\n\t--- Lower Incoming Electron Energy ",energy_0," ---"
    #         print "half energy        =", energy_0/2.0
    #         print "max knock energy   =", (energy_0-binding_energy)/2.0
    #         print "max sampled energy =", recoil_0[len(recoil_0) -1]
    #         print "nudged max knock   =", (energy_0-binding_energy)/2.0*(1.0+1e-5)
    #         print "difference in max sampled = ", (energy_0-binding_energy)/2.0 - recoil_0[len(recoil_0) -1]
    #         outgoing_E_0, angle_0 = linlinlin_dist.sample( energy_0 )
    #         linlinlog_outgoing_E_0, linlinlog_angle_0 = linlinlog_dist.sample( energy_0 )
    #         logloglog_outgoing_E_0, logloglog_angle_0 = logloglog_dist.sample( energy_0 )
    #         print "\tLin-Lin-Lin: E_0 = ",'%.18e' % outgoing_E_0,"\tangle_0 = ",'%.18e' % angle_0
    #         print "\tLin-Lin-Log: E_0 = ",'%.18e' % linlinlog_outgoing_E_0,"\tangle_0 = ",'%.18e' % linlinlog_angle_0

    #         print "\tLog-Log-Log: E_0 = ",'%.18e' % logloglog_outgoing_E_0,"\tangle_0 = ",'%.18e' % logloglog_angle_0

    #         print "\n\t--- Upper Incoming Electron Energy ",energy_1," ---"
    #         print "half energy      =", energy_1/2.0
    #         print "max knock energy =", (energy_1-binding_energy)/2.0
    #         print "max sampled energy =", recoil_1[len(recoil_1) -1]
    #         print "nudged max knock   =", (energy_1-binding_energy)/2.0*(1.0+1e-5)
    #         print "difference in max sampled = ", (energy_1-binding_energy)/2.0 - recoil_1[len(recoil_1) -1]
    #         outgoing_E_1, angle_1 = linlinlin_dist.sample( energy_1 )
    #         linlinlog_outgoing_E_1, linlinlog_angle_1 = linlinlog_dist.sample( energy_1 )
    #         logloglog_outgoing_E_1, logloglog_angle_1 = logloglog_dist.sample( energy_1 )
    #         print "\tLin-Lin-Lin: E_1 = ",'%.18e' % outgoing_E_1,"\tangle_1 = ",'%.18e' % angle_1
    #         print "\tLin-Lin-Log: E_1 = ",'%.18e' % linlinlog_outgoing_E_1,"\tangle_1 = ",'%.18e' % linlinlog_angle_1

    #         print "\tLog-Log-Log: E_1 = ",'%.18e' % logloglog_outgoing_E_1,"\tangle_1 = ",'%.18e' % logloglog_angle_1

    #         logloglog_outgoing_E_1, logloglog_angle_1 = logloglog_dist.sample( energy_1 )

    #         E_in_log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
    #         E_in_lin_interp = (energy-energy_0)/(energy_1-energy_0)

    #         loglog_outgoing_E = numpy.exp(numpy.log(outgoing_E_0) + numpy.log(outgoing_E_1/outgoing_E_0)*E_in_log_interp)
    #         linlog_outgoing_E = outgoing_E_0 + (outgoing_E_1 - outgoing_E_0)*E_in_log_interp
    #         linlin_outgoing_E = outgoing_E_0 + (outgoing_E_1 - outgoing_E_0)*E_in_lin_interp

    #         print "\n\t--- Incoming Electron Energy ",energy," ---"
    #         print "\tLinLin Interp Outgoing Energy: ",'%.18e' % linlin_outgoing_E
    #         print "\tLinLog Interp Outgoing Energy: ",'%.18e' % linlog_outgoing_E
    #         print "\tLogLog Interp Outgoing Energy: ",'%.18e' % loglog_outgoing_E

    #         sample_linlinlin_E, sample_linlinlin_angle = linlinlin_dist.sample( energy )
    #         sample_linlinlog_E, sample_linlinlog_angle = linlinlog_dist.sample( energy )
    #         sample_logloglog_E, sample_logloglog_angle = logloglog_dist.sample( energy )
    #         print "\tLin-Lin-Lin: Outgoing Energy: = ",'%.18e' % sample_linlinlin_E,"\tangle = ",'%.18e' % sample_linlinlin_angle
    #         print "\tLin-Lin-Log: Outgoing Energy: = ",'%.18e' % sample_linlinlog_E,"\tangle = ",'%.18e' % sample_linlinlog_angle
    #         print "\tLog-Log-Log: Outgoing Energy: = ",'%.18e' % sample_logloglog_E,"\tangle = ",'%.18e' % sample_logloglog_angle


    #         energies = [energy_0, energy, energy_1, 2.0*1.49753752578145799e-05 + binding_energy]
    #         for energy in energies:
    #             sample_linlinlin_E_out, sample_linlinlin_E_knock, sample_linlinlin_angle_out, sample_linlinlin_angle_knock = linlinlin_dist.samplePrimaryAndSecondary( energy )
    #             sample_linlinlog_E_out, sample_linlinlog_E_knock, sample_linlinlog_angle_out, sample_linlinlog_angle_knock = linlinlog_dist.samplePrimaryAndSecondary( energy )
    #             sample_logloglog_E_out, sample_logloglog_E_knock, sample_logloglog_angle_out, sample_logloglog_angle_knock = logloglog_dist.samplePrimaryAndSecondary( energy )