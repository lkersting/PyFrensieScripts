#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Data.ACE as ACE
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

# -------------------------------------------------------------------------- ##
#  Electroionization Data
# -------------------------------------------------------------------------- ##
h_native_list = cs_list.get( 'H-Native' )
h_ace_list = cs_list.get( 'H_v14' )
file_names = ['/home/lkersting/frensie/src/packages/test_files/ace/test_pb_epr_ace_file.txt', '/home/lkersting/frensie/src/packages/test_files/ace/test_pb_epr14_ace_file.txt']
table_names = ['82000.12p','82000.14p']

### -------------------------------------------------------------------------- ##
###  Forward Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
h_native_file_name = datadir + h_native_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( h_native_file_name )
h_ace_file_name = datadir + h_ace_list.get( 'electroatomic_file_path' )
ace_file = ACE.ACEFileHandler( h_ace_file_name, "1000.14p", 1 )
ace_data = ACE.XSSEPRDataExtractor( ace_file.getTableNXSArray(), ace_file.getTableJXSArray(), ace_file.getTableXSSArray() )

energy_grid = native_data.getElectronEnergyGrid()

native_dist = Collision.createScreenedRutherfordElasticDistribution( native_data )
ace_dist = Collision.createScreenedRutherfordElasticDistribution( ace_data )

energies = [1e5, 1e-4]
angles = [0.999999, 0.9999995, 1.0]

print "\n----------------------------"
print "----- Native H Tests -----"
print "----------------------------"

print "\n--- evaluate ---\n";

for energy in energies:
  print "energy =", energy
  for angle in angles:
    cdf = native_dist.evaluate( energy, angle )
    print "\tevaluate[",'%.6e' % energy,",",'%.10e' % angle,"] =",'%.16e' % cdf

print "\n--- evaluatePDF ---\n";

for energy in energies:
  print "energy =", energy
  for angle in angles:
    cdf = native_dist.evaluatePDF( energy, angle )
    print "\tevaluatePDF[",'%.6e' % energy,",",'%.10e' % angle,"] =",'%.16e' % cdf

print "\n--- evaluateCDF ---\n";

for energy in energies:
  print "energy =", energy
  for angle in angles:
    cdf = native_dist.evaluateCDF( energy, angle )
    print "\tevaluateCDF[",'%.6e' % energy,",",'%.10e' % angle,"] =",'%.16e' % cdf


print "\n--- sample ---\n";

random_numbers = [0.5]

for energy in energies:
  Prng.RandomNumberGenerator.setFakeStream(random_numbers)
  e_out, scattering_angle_cosine = native_dist.sample( energy )
  print "\tenergy =",energy,"\trandom_number =",random_numbers[0]
  print "scattering mu =",'%.16e' % scattering_angle_cosine, "\toutgoing energy =",'%.16e' % e_out

print "\n----------------------------"
print "----- Ace H Tests -----"
print "----------------------------"

print "\n--- evaluate ---\n";

for energy in energies:
  print "energy =", energy
  for angle in angles:
    cdf = ace_dist.evaluate( energy, angle )
    print "\tevaluate[",'%.6e' % energy,",",'%.10e' % angle,"] =",'%.16e' % cdf

print "\n--- evaluatePDF ---\n";

for energy in energies:
  print "energy =", energy
  for angle in angles:
    cdf = ace_dist.evaluatePDF( energy, angle )
    print "\tevaluatePDF[",'%.6e' % energy,",",'%.10e' % angle,"] =",'%.16e' % cdf

print "\n--- evaluateCDF ---\n";

for energy in energies:
  print "energy =", energy
  for angle in angles:
    cdf = ace_dist.evaluateCDF( energy, angle )
    print "\tevaluateCDF[",'%.6e' % energy,",",'%.10e' % angle,"] =",'%.16e' % cdf


print "\n--- sample ---\n";

random_numbers = [0.5]

for energy in energies:
  Prng.RandomNumberGenerator.setFakeStream(random_numbers)
  e_out, scattering_angle_cosine = ace_dist.sample( energy )
  print "\tenergy =",energy,"\trandom_number =",random_numbers[0]
  print "scattering mu =",'%.16e' % scattering_angle_cosine, "\toutgoing energy =",'%.16e' % e_out