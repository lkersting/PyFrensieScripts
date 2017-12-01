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

# -------------------------------------------------------------------------- ##
#  Adjoint Elastic Data
# -------------------------------------------------------------------------- ##
elements = ['H-Native']
interps = ["LogLogLog", "LinLinLin", "LinLinLog"]
interps = ["LogLogLog"]
energies = [1e-5, 1e-3, 20.0 ]

for z in elements:
    print "\n----------------------------"
    print "-----", z, "Adjoint Tests -----"
    print "----------------------------"
    data_list = cs_list.get( z )

    adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
    adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
    energy_grid = adjoint_data.getAdjointElectronEnergyGrid()


    tot_adjoint_elastic_cs = adjoint_data.getAdjointTotalElasticCrossSection()
    adjoint_cutoff_cs = adjoint_data.getAdjointCutoffElasticCrossSection()
    reduced_cutoff_ratio = adjoint_data.getReducedCutoffCrossSectionRatios()
    adjoint_screen_rutherford_cs = adjoint_data.getAdjointScreenedRutherfordElasticCrossSection()
    adjoint_screen_rutherford_index = adjoint_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
    moment_cs = adjoint_data.getAdjointMomentPreservingCrossSection()
    adjoint_moment_index = adjoint_data.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex()

    ###
    ###  Moment Preserving Adjoint Distribution/Reaction Unit Test Check
    ###
    for interp in interps:
        print "\n--- ",interp,"Tests ---"

        for energy in energies:

            index = 0
            for i in range(0, energy_grid.size ):
                if energy_grid[i] <= energy:
                    index = i

            energy_0 = energy_grid[index]
            cs_0 = moment_cs[index]
            cs = cs_0

            if energy_0 != energy:
                energy_1 = energy_grid[index+1]
                cs_1 = moment_cs[index+1]

                log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
                lin_interp = (energy-energy_0)/(energy_1-energy_0)

                cs = cs_0 + (cs_1-cs_0)*lin_interp

            print '\tEnergy = ', '%.6e' % energy,'\tcs  = ','%.16e' % cs


