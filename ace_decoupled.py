#! /usr/bin/env python
import PyFrensie.Data as Data
import PyFrensie.Data.ACE as ACE
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.MonteCarlo as MonteCarlo
import PyFrensie.MonteCarlo.Electron as Electron
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

Utility.initFrensiePrng()

datadir = '/home/software/mcnpdata/'
database_path = datadir + 'database.xml'

database = Data.ScatteringCenterPropertiesDatabase(database_path)
au_properties = database.getAtomProperties( Data.ZAID(79000) )

au_electron_prop = au_properties.getSharedElectroatomicDataProperties(
                            Data.ElectroatomicDataProperties.ACE_EPR_FILE, 14 )


file_name = datadir + au_electron_prop.filePath()
table_start = au_electron_prop.fileStartLine()
table_name = au_electron_prop.tableName()


### -------------------------------------------------------------------------- ##
###  ACE Elastic Check
### -------------------------------------------------------------------------- ##
ace_file = ACE.ACEFileHandler( file_name, table_name, table_start )
xss_extractor = ACE.XSSEPRDataExtractor( ace_file.getTableNXSArray(), ace_file.getTableJXSArray(), ace_file.getTableXSSArray() )

decoupled_react = Electron.createDecoupledElasticReaction( xss_extractor )


###
###  Decoupled Distribution/Reaction Check
###

energies = [15.7, 10.0]

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
  if sampling_ratio > 0.5 and sampling_ratio < 1.0:
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
    shell = decoupled_react.react( electron, bank )
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
    shell = decoupled_react.react( electron, bank )
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
    print '\tangle[','%.6e' %cdfs[i],']\t= ','%.16e' % angles[i]
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
    shell = decoupled_react.react( electron, bank )
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