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
#  Elastic Data
# -------------------------------------------------------------------------- ##
elements = ['Pb-Native','Al-Native','H-Native']
interps = ["LogLogLog", "LinLinLin", "LinLinLog"]
energies = [1e-5, 1e-3, 1e5 ]

for z in elements:
    print "\n----------------------------"
    print "-----", z, "Tests -----"
    print "----------------------------"
    data_list = cs_list.get( z )
    file_name = datadir + data_list.get( 'electroatomic_file_path' )
    native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
    energy_grid = native_data.getElectronEnergyGrid()


    tot_elastic_cs = native_data.getTotalElasticCrossSection()
    cutoff_cs = native_data.getCutoffElasticCrossSection()
    screen_rutherford_cs = native_data.getScreenedRutherfordElasticCrossSection()
    screen_rutherford_index = native_data.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
    moment_cs = native_data.getMomentPreservingCrossSection()
    moment_index = native_data.getMomentPreservingCrossSectionThresholdEnergyIndex()

    print screen_rutherford_index
    print screen_rutherford_cs.size
    print energy_grid.size
    if z == 'Pb-Native':
        energies = [1e-5,4e-4,1e5]
    else:
        energies = [1e-5, 1e-3, 1e5, 4.3750e1 ]

    for interp in interps:
        print "\n--- ",interp,"Tests ---"
        print "\t--- Moment Preserving Tests ---"

        ###
        ###  Moment Preserving Reaction Unit Test Check
        ###
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

                cs = 0.0
                if interp =="LogLogLog":
                    cs = (cs_0)*pow((cs_1)/(cs_0),log_interp)
                elif interp == "LinLinLin":
                    cs = cs_0 + (cs_1-cs_0)*lin_interp
                else:
                    cs = cs_0 + (cs_1-cs_0)*log_interp

            print '\tEnergy = ',energy,'\tcs  = ','%.16e' % cs

        print "\t--- SR Elastic Tests ---"
        ###
        ###  SR Elastic Reaction Unit Test Check
        ###
        for energy in energies:

            index = 0
            for i in range(0, energy_grid.size ):
                if energy_grid[i] <= energy:
                    index = i

            energy_0 = energy_grid[index]
            tot_cs_0 = tot_elastic_cs[index]
            cut_cs_0 = cutoff_cs[index]
            cs_0 = 0
            sr_index = index - screen_rutherford_index
            if sr_index >= 0:
                cs_0 = screen_rutherford_cs[sr_index]

            cs = cs_0
            tot_cs = tot_cs_0
            cut_cs = cut_cs_0

            if energy_0 != energy:
                energy_1 = energy_grid[index+1]
                cs_1 = screen_rutherford_cs[sr_index+1]
                tot_cs_1 = tot_elastic_cs[index+1]
                cut_cs_1 = cutoff_cs[index+1]

                log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
                lin_interp = (energy-energy_0)/(energy_1-energy_0)

                cs = 0.0
                tot_cs = 0.0
                cut_cs = 0.0
                if interp =="LogLogLog":
                    cs = (cs_0)*pow((cs_1)/(cs_0),log_interp)
                    tot_cs = (tot_cs_0)*pow((tot_cs_1)/(tot_cs_0),log_interp)
                    cut_cs = (cut_cs_0)*pow((cut_cs_1)/(cut_cs_0),log_interp)
                elif interp == "LinLinLin":
                    cs = cs_0 + (cs_1-cs_0)*lin_interp
                    tot_cs = tot_cs_0 + (tot_cs_1-tot_cs_0)*lin_interp
                    cut_cs = cut_cs_0 + (cut_cs_1-cut_cs_0)*lin_interp
                else:
                    cs = cs_0 + (cs_1-cs_0)*log_interp
                    tot_cs = tot_cs_0 + (tot_cs_1-tot_cs_0)*log_interp
                    cut_cs = cut_cs_0 + (cut_cs_1-cut_cs_0)*log_interp

            sr_cs = tot_cs - cut_cs
            print '\tEnergy = ',energy,'\tcs (interp) = ','%.16e' % cs
            print '\tEnergy = ',energy,'\tcs (calc)   = ','%.16e' % sr_cs