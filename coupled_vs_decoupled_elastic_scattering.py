#! /usr/bin/env python
import PyFrensie.Data as Data
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.MonteCarlo as MonteCarlo
import PyFrensie.MonteCarlo.Electron as Electron
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

Utility.initFrensiePrng()

#datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/frensie/src/packages/test_files/'
name = 'native/test_epr_1_native.xml'

# database = datadir + '/cross_sections.xml'
# database = Data.ScatteringCenterPropertiesDatabase(database)
# h_properties = database.getAtomProperties( Utility.ZAID(1001) )
# data_list = cs_list.get( 'H-Native' )

### -------------------------------------------------------------------------- ##
###  Forward Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
native_file_name = datadir + name
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
methods = [Electron.SIMPLIFIED_UNION,Electron.ONE_D_UNION,Electron.TWO_D_UNION]
energies = [1.0e+5,6.625E+01,200.0]

interps = ["LogLogLog"]
methods = [Electron.SIMPLIFIED_UNION]
energies = [0.314, 0.521]
for interp in interps:
  print "\n-----", interp ,"-----\n"

  coupled_react = None
  coupled_dist = None
  decoupled_react = None
  cutoff_dist = None
  for method in methods:
      if interp == "LinLinLog":
        # Coupled Reaction
        coupled_react = Electron.createCoupledElasticReaction_LinLogCorrelated( native_data, method, 1e-15 )
        # Coupled Distribution
        coupled_dist = Electron.createCoupledElasticDistribution_LinLogCorrelated( native_data, method, 1e-15 )
        # Decoupled Reaction
        decoupled_react = Electron.createDecoupledElasticReaction_LinLogCorrelated( native_data, 1e-15 )
        # Cutoff Distribution
        cutoff_dist = Electron.createCutoffElasticDistribution_LinLogCorrelated( native_data, 0.999999, 1e-15 )
      if interp == "LogLogLog":
        # Coupled Reaction
        coupled_react = Electron.createCoupledElasticReaction_LogLogCorrelated( native_data, method, 1e-15 )
        # Coupled Distribution
        coupled_dist = Electron.createCoupledElasticDistribution_LogLogCorrelated( native_data, method, 1e-15 )
        # Decoupled Reaction
        decoupled_react = Electron.createDecoupledElasticReaction_LogLogCorrelated( native_data, 1e-15 )
        # Cutoff Distribution
        cutoff_dist = Electron.createCutoffElasticDistribution_LogLogCorrelated( native_data, 0.999999, 1e-15 )
      elif interp == "LinLinLin":
        # Coupled Reaction
        coupled_react = Electron.createCoupledElasticReaction_LinLinCorrelated( native_data, method, 1e-15 )
        # Coupled Distribution
        coupled_dist = Electron.createCoupledElasticDistribution_LinLinCorrelated( native_data, method, 1e-15 )
        # Decoupled Reaction
        decoupled_react = Electron.createDecoupledElasticReaction_LinLinCorrelated( native_data, 1e-15 )
        # Cutoff Distribution
        cutoff_dist = Electron.createCutoffElasticDistribution_LinLinCorrelated( native_data, 0.999999, 1e-15 )


      print "\n----- ",method," -----\n"

      n = 100
      electron = MonteCarlo.ElectronState( 0 )
      bank = MonteCarlo.ParticleBank()

      for energy in energies:
        fig, ax = plt.subplots()

        print "\t -- React 256 keV--"
        energy0 = 0.256
        electron.setEnergy( energy0 )
        sampling_ratio = decoupled_react.getSamplingRatio(energy0)
        print "sampling_ratio = ", sampling_ratio
        random_numbers = [0.0]*2*n
        if sampling_ratio > 0.5:
          cdf1 = numpy.linspace(0.5, sampling_ratio-1e-15, n/2)
          cdf2 = numpy.linspace(sampling_ratio+1e-15, 1.0-1e-15, n/2)
          cdfs = numpy.append(cdf1, cdf2 )
        else:
          cdfs = numpy.linspace(0.5, 1.0-1e-15, n)

        for i in range(n):
          j = 2*i
          random_numbers[j] = cdfs[i]
        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

        angles = [None]*n
        for i in range(n):
          # print random_numbers[2*i], random_numbers[2*i+1]
          electron.setDirection( 0.0, 0.0, 1.0 )
          shell = coupled_react.react( electron, bank )
          angles[i] = electron.getZDirection()
          # print '\tangle[','%.6e' %cdfs[i],']\t= ','%.16e' % angles[i]
        Prng.RandomNumberGenerator.unsetFakeStream()
        ax.semilogy(cdfs, angles, label="256 keV")


        print "Energy = ",energy
        electron = MonteCarlo.ElectronState( 0 )
        electron.setEnergy( energy )

        print "\t -- React Coupled--"
        sampling_ratio = decoupled_react.getSamplingRatio(energy)
        print "sampling_ratio = ", sampling_ratio
        random_numbers = [0.0]*2*n
        if sampling_ratio > 0.5:
          cdf1 = numpy.linspace(0.5, sampling_ratio-1e-15, n/2)
          cdf2 = numpy.linspace(sampling_ratio+1e-15, 1.0-1e-15, n/2)
          cdfs = numpy.append(cdf1, cdf2 )
        else:
          cdfs = numpy.linspace(0.5, 1.0-1e-15, n)
        for i in range(n):
          j = 2*i
          random_numbers[j] = cdfs[i]
        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

        angles = [None]*n
        for i in range(n):
          # print random_numbers[2*i], random_numbers[2*i+1]
          electron.setDirection( 0.0, 0.0, 1.0 )
          shell = coupled_react.react( electron, bank )
          angles[i] = electron.getZDirection()
          # print '\tangle[','%.6e' %cdfs[i],']\t= ','%.16e' % angles[i]
        Prng.RandomNumberGenerator.unsetFakeStream()
        ax.semilogy(cdfs, angles, label=str(int(energy*1000)) +" keV Coupled")

        print "\t -- React Decoupled--"
        sampling_ratio = decoupled_react.getSamplingRatio(energy)
        print "sampling_ratio = ", sampling_ratio
        random_numbers = [0.0]*3*n
        for i in range(n):
          # Set the random number to sample between cutoff or peak
          j = 3*i
          random_numbers[j] = 0.0#cdfs[i]
          # Set the random number to sample the distribution
          k = j+1
          if cdfs[i] <= sampling_ratio:
            random_numbers[j] = 0.0#cdfs[i]
            random_numbers[k] = cdfs[i]/sampling_ratio
          else:
            random_numbers[j] = 1.0-1e-15#cdfs[i]
            random_numbers[k] = (cdfs[i] - sampling_ratio)/(1.0 - sampling_ratio)
        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

        angles = [None]*n
        for i in range(n):
          # print random_numbers[3*i], random_numbers[3*i+1], random_numbers[3*i+2]
          electron.setDirection( 0.0, 0.0, 1.0 )
          shell = decoupled_react.react( electron, bank )
          angles[i] = electron.getZDirection()
          # print '\tangle[','%.6e' %cdfs[i],']\t= ','%.16e' % angles[i]
        ax.semilogy(cdfs, angles, label=str(int(energy*1000)) +" keV Decoupled")

        print "\t -- React 10 MeV Decoupled--"
        energy1 = 10.0
        electron.setEnergy( energy1 )
        sampling_ratio = decoupled_react.getSamplingRatio(energy1)
        print "sampling_ratio = ", sampling_ratio
        random_numbers = [0.0]*3*n
        for i in range(n):
          # Set the random number to sample between cutoff or peak
          j = 3*i
          random_numbers[j] = 0.0
          # Set the random number to sample the distribution
          k = j+1
          random_numbers[k] = cdfs[i]

        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

        angles = [None]*n
        for i in range(n):
          # print random_numbers[3*i], random_numbers[3*i+1], random_numbers[3*i+2]
          electron.setDirection( 0.0, 0.0, 1.0 )
          shell = decoupled_react.react( electron, bank )
          angles[i] = electron.getZDirection()
          # print '\tangle[','%.6e' %cdfs[i],']\t= ','%.16e' % angles[i]
        ax.semilogy(cdfs, angles, label="10 MeV Cutoff")

        print "\t -- React 10 MeV Coupled--"
        energy = 10.0
        random_numbers = [0.0]*2*n
        if sampling_ratio > 0.5:
          cdf1 = numpy.linspace(0.5, sampling_ratio-1e-15, n/2)
          cdf2 = numpy.linspace(sampling_ratio+1e-15, 1.0-1e-15, n/2)
          cdfs = numpy.append(cdf1, cdf2 )
        else:
          cdfs = numpy.linspace(0.5, 1.0-1e-15, n)

        for i in range(n):
          j = 2*i
          random_numbers[j] = cdfs[i]
        Prng.RandomNumberGenerator.setFakeStream(random_numbers)

        angles = [None]*n
        for i in range(n):
          # print random_numbers[2*i], random_numbers[2*i+1]
          electron.setDirection( 0.0, 0.0, 1.0 )
          shell = coupled_react.react( electron, bank )
          angles[i] = electron.getZDirection()
          # print '\tangle[','%.6e' %cdfs[i],']\t= ','%.16e' % angles[i]
        Prng.RandomNumberGenerator.unsetFakeStream()
        ax.semilogy(cdfs, angles, label="10 MeV Coupled")


        plt.xlim(0.5, 1.0)
        plt.ylim(0.99995, 1.000001)
        plt.yticks(numpy.arange(0.99995, 1.000001, 0.00001))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.5f'))
        plt.xlabel("CDF")
        plt.ylabel("Sampled Angle Cosine")
        plt.title("Differences between Coupled and Decoupled Elastic")

        plt.legend(loc=4)
        plt.show()