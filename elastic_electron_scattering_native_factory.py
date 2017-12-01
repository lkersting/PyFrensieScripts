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

### -------------------------------------------------------------------------- ##
###  Forward Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
name = 'Pb-Native'
data_list = cs_list.get( name )
native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )

interps = ["LinLinLin", "LinLinLog"]
energies = [1e-4]

for interp in interps:
    print "\n----------------------------"
    print "--- ",interp,name,"Tests ---"
    print "----------------------------"

    hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, interp, True, 1e-14)

    ###
    ###  Sample
    ###
    print '\n--- Sampling ---'

    for energy in energies:
        print "\n\t Energy = ", energy
        random_numbers = [0.5, 0.9, 0.95, 1.0-1e-15]
        Prng.RandomNumberGenerator.setFakeStream(random_numbers)
        for rnd in random_numbers:
            energy, angle = hybrid_dist.sample( energy )
            print "\tsample[",energy,",",'%.16e' % rnd,"] = ",'%.16e' % angle


### -------------------------------------------------------------------------- ##
###  Adjoint Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
name = 'H-Native'
data_list = cs_list.get( name )
native_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
native_data = Native.AdjointElectronPhotonRelaxationDataContainer( native_file_name )
energy_grid = native_data.getAdjointElectronEnergyGrid()

tot_elastic_cs = native_data.getAdjointTotalElasticCrossSection()
cutoff_cs = native_data.getAdjointCutoffElasticCrossSection()
screen_rutherford_cs = native_data.getAdjointScreenedRutherfordElasticCrossSection()
screen_rutherford_index = native_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
moment_cs = native_data.getAdjointMomentPreservingCrossSection()
moment_index = native_data.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex()
discrete_energy_grid = native_data.getAdjointElasticAngularEnergyGrid()

interps = ["LinLinLog"]
energies = [1e-3]

for interp in interps:
    print "\n----------------------------"
    print "--- ",interp,name,"Adjoint Tests ---"
    print "----------------------------"

    hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, interp, True, 1e-14)

    ###
    ###  Sample
    ###
    print '\n--- Sampling ---'

    for energy in energies:
        print "\n\t Energy = ", energy
        random_numbers = [0.1, 0.2, 0.4, 1.0-1e-15]
        Prng.RandomNumberGenerator.setFakeStream(random_numbers)
        for rnd in random_numbers:
            energy, angle = hybrid_dist.sample( energy )
            print "\tsample[",energy,",",'%.16e' % rnd,"] = ",'%.16e' % angle

