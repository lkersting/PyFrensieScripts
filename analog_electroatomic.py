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
data_list = cs_list.get( 'Pb-Native' )

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
###  Lin-Log Analog Reaction Unit Test Check
###
print "\n----- Lin-Lin-Log -----"
analog_dist = Collision.createAnalogElasticDistribution(native_data, "LinLinLog", True, 1e-15)

energies = [1e-5,4e-4,1e5]
for energy in energies:
    print "Energy = ",energy

    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    energy_0 = energy_grid[index]
    cs_0 = cutoff_cs[index]
    cs = cs_0

    if energy_0 != energy:
        energy_1 = energy_grid[index+1]
        cs_1 = cutoff_cs[index+1]
        cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )

    print '\tcutoff cs  = ','%.16e' % cs

    cutoff_cdf = analog_dist.evaluateCDFAtCutoff( energy )
    print '\tcutoff_cdf = ','%.16e' % cutoff_cdf

    total_cs = cs/cutoff_cdf
    print '\ttotal cs   = ','%.16e' % total_cs

###
###  Lin-Lin Analog Reaction Unit Test Check
###
print "\n----- Lin-Lin-Lin -----"
analog_dist = Collision.createAnalogElasticDistribution(native_data, "LinLinLin", True, 1e-15)

energies = [1e-5,4e-4,1e5]
for energy in energies:
    print "Energy = ",energy

    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    energy_0 = energy_grid[index]
    cs_0 = cutoff_cs[index]
    cs = cs_0

    if energy_0 != energy:
        energy_1 = energy_grid[index+1]
        cs_1 = cutoff_cs[index+1]
        cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )

    print '\tcutoff cs  = ','%.16e' % cs

    cutoff_cdf = analog_dist.evaluateCDFAtCutoff( energy )
    print '\tcutoff_cdf = ','%.16e' % cutoff_cdf

    total_cs = cs/cutoff_cdf
    print '\ttotal cs   = ','%.16e' % total_cs

