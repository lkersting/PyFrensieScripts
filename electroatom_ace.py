#! /usr/bin/env python
import PyFrensie.Data as Data
import PyFrensie.Data.ACE as ACE
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.MonteCarlo.Collision as Collision
import numpy
import matplotlib.pyplot as plt

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


# -------------------------------------------------------------------------- ##
# Electroatom Tests
# -------------------------------------------------------------------------- ##

ace_file = ACE.ACEFileHandler( file_name, table_name, table_start )
xss_extractor = ACE.XSSEPRDataExtractor( ace_file.getTableNXSArray(), ace_file.getTableJXSArray(), ace_file.getTableXSSArray() )

energy_grid = xss_extractor.extractElectronEnergyGrid()


###
### Electroatom/Electroatom Core Test Check
###
print "\n----- Electroatom ACE Factory Class -----"
brem_cross_sections = xss_extractor.extractBremsstrahlungCrossSection()
excitation_cross_sections = xss_extractor.extractExcitationCrossSection()
ionization_cross_sections = xss_extractor.extractElectroionizationCrossSection()
tot_elastic_cross_sections = xss_extractor.extractElasticTotalCrossSection()
cutoff_cross_sections = xss_extractor.extractElasticCutoffCrossSection()

energies = [2e-3, 4e-4, 9e-5, 1e-5, 1e5]
for energy in energies:

    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    print "\nEnergy = ",energy,'\tindex = ', index

    brem_cs = brem_cross_sections[index]
    excitation_cs = excitation_cross_sections[index]
    cutoff_cs = cutoff_cross_sections[index]
    analog_cs = tot_elastic_cross_sections[index]
    ionization_cs = ionization_cross_sections[index]


    inelastic_cs = brem_cs + excitation_cs + ionization_cs
    total_cs_analog = inelastic_cs + analog_cs
    total_cs_cutoff = inelastic_cs + cutoff_cs

    print '\tbrem_cs       = ','%.16e' % brem_cs
    print '\texcitation_cs = ','%.16e' % excitation_cs
    print '\tionization_cs = ','%.16e' % ionization_cs
    print '\t------------------------------------------------'
    print '\tinelastic_cs  = ','%.16e' % inelastic_cs
    print '\t------------------------------------------------'
    print '\tanalog_cs (lin)  = ','%.16e' % analog_cs
    print '\tcutoff_cs (log) = ','%.16e' % cutoff_cs
    print '\t------------------------------------------------'
    print '\ttotal cs (analog) = ','%.16e' % total_cs_analog
    print '\ttotal cs (cutoff) = ','%.16e' % total_cs_cutoff

energies = [1.0e5, 1.995260e1, 6.309570e0, 1.995260e-3, 1.995260e-4, 1.0e-5]
print "\n--- Analog Cross Section ---"
for energy in energies:
    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    energy_0 = energy_grid[index]
    elastic_cs = 0.0
    if energy_0 != energy:
        energy_1 = energy_grid[index+1]
        lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )

        elastic_cs_0 = cutoff_cross_sections[index]
        elastic_cs_1 = cutoff_cross_sections[index+1]

        elastic_cs = elastic_cs_0 + (elastic_cs_1 - elastic_cs_0)*lin_interp
    else:
        elastic_cs = cutoff_cross_sections[index]
    print '\tcs[','%.6e' %energy,']:','%.16e' % elastic_cs

energies = [1.0e5, 1.995260e1, 6.30957, 1e-3, 1.995260e-4, 1.0e-5]
print "\n--- Cutoff Cross Section ---"

for energy in energies:
    print "\nEnergy = ",energy
    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    cutoff_cs = cutoff_cross_sections[index]
    print '\tcutoff_cs','%.16e' % cutoff_cs


print "\n--- Electroatom Factory ---"
energies = [2e-3, 4e-4, 9e-5, 1e-5, 1e5]
for energy in energies:

   index = 0
   for i in range(0, energy_grid.size ):
       if energy_grid[i] == energy:
           index = i

   print "\nEnergy = ",energy,'\tindex = ', index

   brem_cs = brem_cross_sections[index]
   excitation_cs = excitation_cross_sections[index]
   cutoff_cs = cutoff_cross_sections[index]
   analog_cs = tot_elastic_cross_sections[index]
   ionization_cs = ionization_cross_sections[index]

   inelastic_cs = brem_cs + excitation_cs + ionization_cs
   total_cs = inelastic_cs + analog_cs

   print '\tbrem_cs       = ','%.16e' % brem_cs
   print '\texcitation_cs = ','%.16e' % excitation_cs
   print '\tionization_cs = ','%.16e' % ionization_cs
   print '\t-----------------------------------------'
   print '\tinelastic_cs  = ','%.16e' % inelastic_cs
   print '\telastic_cs    = ','%.16e' % analog_cs
   print '\t-----------------------------------------'
   print '\ttotal cs      = ','%.16e' % total_cs

###
### Electroatom ACE Factory/Electroatom Factory Test Check
###
print "\n----- Electroatom -----\n"
print "\n----- Au -----\n"

ace_file = ACE.ACEFileHandler( file_name, table_name, table_start )
xss_extractor = ACE.XSSEPRDataExtractor( ace_file.getTableNXSArray(), ace_file.getTableJXSArray(), ace_file.getTableXSSArray() )

energy_grid = xss_extractor.extractElectronEnergyGrid()

brem_cross_sections = xss_extractor.extractBremsstrahlungCrossSection()
excitation_cross_sections = xss_extractor.extractExcitationCrossSection()
ionization_cross_sections = xss_extractor.extractElectroionizationCrossSection()
tot_elastic_cross_sections = xss_extractor.extractElasticTotalCrossSection()
cutoff_cross_sections = xss_extractor.extractElasticCutoffCrossSection()

energies = [5.2371421547030929e-02, 2e-3, 4e-4, 9e-5, 1e-5, 1e5]
for energy in energies:
  index = 0
  for i in range(0, energy_grid.size ):
    if energy_grid[i] <= energy:
        index = i
  print energy_grid[index]
  if energy != energy_grid[index]:
    print energy_grid[index+1]

    energy_0 = energy_grid[index]
    brem_cs_0 = brem_cross_sections[index]
    excitation_cs_0 = excitation_cross_sections[index]
    ionization_cs_0 = ionization_cross_sections[index]
    analog_cs_0 = tot_elastic_cross_sections[index]
    cs_0 = brem_cs_0 + excitation_cs_0 + ionization_cs_0 + analog_cs_0

    energy_1 = energy_grid[index+1]
    brem_cs_1 = brem_cross_sections[index+1]
    excitation_cs_1 = excitation_cross_sections[index+1]
    ionization_cs_1 = ionization_cross_sections[index+1]
    analog_cs_1 = tot_elastic_cross_sections[index+1]
    cs_1 = brem_cs_1 + excitation_cs_1 + ionization_cs_1 + analog_cs_1


    lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
    log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

    numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)

    brem_cs = numpy.exp(numpy.log(brem_cs_0) + numpy.log(brem_cs_1/brem_cs_0)*log_interp)

    excitation_cs = numpy.exp(numpy.log(excitation_cs_0) + numpy.log(excitation_cs_1/excitation_cs_0)*log_interp)

    ionization_cs = numpy.exp(numpy.log(ionization_cs_0) + numpy.log(ionization_cs_1/ionization_cs_0)*log_interp)

    analog_cs = numpy.exp(numpy.log(analog_cs_0) + numpy.log(analog_cs_1/analog_cs_0)*log_interp)

    tot_cs = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)
  else:
    brem_cs = brem_cross_sections[index]
    excitation_cs = excitation_cross_sections[index]
    ionization_cs = ionization_cross_sections[index]
    analog_cs = tot_elastic_cross_sections[index]
  tot_cs = brem_cs + excitation_cs + ionization_cs + analog_cs

  max_excitation = (excitation_cs)/tot_cs
  max_brem = (brem_cs+excitation_cs)/tot_cs
  max_ionization = (ionization_cs+brem_cs+excitation_cs)/tot_cs
  max_elastic = (ionization_cs+brem_cs+excitation_cs+analog_cs)/tot_cs


  print "\nenergy        = ",'%.16e' % energy
  print '\tbrem_cs       = ','%.16e' % brem_cs
  print '\texcitation_cs = ','%.16e' % excitation_cs
  print '\tionization_cs = ','%.16e' % ionization_cs
  print '\telastic_cs    = ','%.16e' % analog_cs
  print '\t--------------------------------'
  print '\ttot_cs        = ','%.16e' % tot_cs


  print '\n\tindex = ', index

  print '\n\tmax excitation random number = ','%.16e' % max_excitation
  print '\tmax brem random number = ','%.16e' % max_brem
  print '\tmax ionization random number = ','%.16e' % max_ionization
  print '\tmax elastic random number = ','%.16e' % max_elastic