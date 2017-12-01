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
energies = [1e-5, 1e-3, 20.0 ]

for z in elements:
    print "\n----------------------------"
    print "-----", z, "Adjoint Tests -----"
    print "----------------------------"
    data_list = cs_list.get( z )

    adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
    adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
    energy_grid = adjoint_data.getAdjointElectronEnergyGrid()


    tot_cs = adjoint_data.getAdjointTotalElasticCrossSection()

    ###
    ###  Decoupled Adjoint Distribution/Reaction Unit Test Check
    ###
    for energy in energies:

        index = 0
        for i in range(0, energy_grid.size ):
            if energy_grid[i] <= energy:
                index = i

        energy_0 = energy_grid[index]
        cs_0 = tot_cs[index]
        cs = cs_0

        if energy_0 != energy:
            energy_1 = energy_grid[index+1]
            cs_1 = tot_cs[index+1]

            log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
            lin_interp = (energy-energy_0)/(energy_1-energy_0)

            cs = cs_0 + (cs_1-cs_0)*lin_interp

        print '\tEnergy = ', '%.6e' % energy,'\tcs  = ','%.16e' % cs


