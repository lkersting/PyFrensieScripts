#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.Utility.Distribution as Distribution
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
import numpy
import matplotlib.pyplot as plt
from collections import Counter

Utility.initFrensiePrng()

#datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/frensie/src/packages/test_files/'

source = Teuchos.FileInputSource( datadir + '/cross_sections.xml' )
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList( xml_obj )

# -------------------------------------------------------------------------- ##
#  Brem Data
# -------------------------------------------------------------------------- ##
# Possible Elements ['H-Native', 'Al-Native', 'Pb-Native']
elements = ['Pb-Native']
# Possible Interpolation Policies ["LogLogLog", "LinLinLin", "LinLinLog"]
interps = ["LogLogLog"]
# Possible Interpolation Schemes ["Unit-base", "Unit-base CDF", "Correlated Unit-base", "Corresponding Energies", "Cumulative Points"]
schemes = ["Unit-base", "Correlated Unit-base", "Corresponding Energies", "Cumulative Points"]
# Show Relative difference in schemes (True/False)
show_difference = False
# Possible energies [1e-2, 1e-1, 1.0, 15.7, 20.0]
energies = [1e3]
# # Step length between plot points
step = 0.0001
length = int(1.0/step)

plot_number = 1
for z in elements:
    print "\n----------------------------"
    print "-----", z, "Tests -----"
    print "----------------------------"
    data_list = cs_list.get(z)
    file_name = datadir + data_list.get('electroatomic_file_path')
    native_data = Native.ElectronPhotonRelaxationDataContainer(file_name)
    energy_grid = native_data.getBremsstrahlungEnergyGrid()

    # Plot the given energy
    for energy in energies:
        print "Energy = ", energy
        print "----------------------------"

        # Find energy in energy grid (lower bin index)
        index = 0
        for i in range(0, len(energy_grid)):
            if energy_grid[i] < energy:
                index = i

        # Get lower and upper energy and the interpolation alpha
        energy_0 = energy_grid[index]
        energy_1 = energy_grid[index+1]
        lin_E_alpha = (energy - energy_0)/(energy_1 - energy_0)
        log_E_alpha = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

        # Get lower energy bin data
        e_losses_0 = native_data.getBremsstrahlungPhotonEnergy(energy_0)
        pdfs_0 = native_data.getBremsstrahlungPhotonPDF(energy_0)
        dist_0 = Distribution.TabularDistribution_LinLin(e_losses_0, pdfs_0)
        log_dist_0 = Distribution.TabularDistribution_LogLog(e_losses_0, pdfs_0)

        # Process lower energy bin data
        log_dist_0 = Distribution.TabularDistribution_LogLog(e_losses_0, pdfs_0)
        log_e_losses_0 = numpy.zeros(shape=(len(pdfs_0)))
        for i in range(0, len(log_e_losses_0)):
            log_e_losses_0[i] = numpy.log(e_losses_0[i])
        process_dist_0 = Distribution.TabularDistribution_LogLin(log_e_losses_0, pdfs_0)

        # Get the lower cdf values
        cdfs_0 = numpy.zeros(shape=(len(pdfs_0)))
        log_cdfs_0 = numpy.zeros(shape=(len(pdfs_0)))
        process_cdfs_0 = numpy.zeros(shape=(len(pdfs_0)))
        for i in range(0, len(log_cdfs_0)):
            log_cdfs_0[i] = log_dist_0.evaluateCDF(e_losses_0[i])
            cdfs_0[i] = dist_0.evaluateCDF(e_losses_0[i])
            diff = log_cdfs_0[i] - cdfs_0[i]
            if log_e_losses_0[i] == numpy.log(e_losses_0[i]):
                print "Diff = (", log_cdfs_0[i], "-", process_cdfs_0[i], ") = ", diff
        for i in range(0, len(log_cdfs_0)):
            log_cdfs_0[i] = log_dist_0.evaluateCDF(e_losses_0[i])
            process_cdfs_0[i] = process_dist_0.evaluateCDF(log_e_losses_0[i])
            diff = log_cdfs_0[i] - process_cdfs_0[i]
            if log_e_losses_0[i] == numpy.log(e_losses_0[i]):
                print "Diff = (", log_cdfs_0[i], "-", process_cdfs_0[i], ") = ", diff
        process_cdf_dist_0 = Distribution.TabularCDFDistribution_LogLin(log_e_losses_0, log_cdfs_0, True)
        for i in range(0, len(log_cdfs_0)):
            log_cdfs_0[i] = log_dist_0.evaluateCDF(e_losses_0[i])
            process_cdfs_0[i] = process_cdf_dist_0.evaluateCDF(log_e_losses_0[i])
            diff = log_cdfs_0[i] - process_cdfs_0[i]
            if log_e_losses_0[i] == numpy.log(e_losses_0[i]):
                print "Diff = (", log_cdfs_0[i], "-", process_cdfs_0[i], ") = ", diff

        e_loss_min = e_losses_0[0]
        e_losses = numpy.logspace(numpy.log10(e_loss_min), numpy.log10(energy_0), num=length)
        log_e_losses = numpy.zeros(shape=(len(e_losses)))
        for i in range(0, len(log_e_losses_0)):
            log_e_losses[i] = numpy.log(e_losses[i])
        for i in range(0, len(log_cdfs_0)):
            log_cdfs_0[i] = log_dist_0.evaluateCDF(e_losses[i])
            process_cdfs_0[i] = process_cdf_dist_0.evaluateCDF(log_e_losses[i])
            diff = log_cdfs_0[i] - process_cdfs_0[i]
            if log_e_losses[i] == numpy.log(e_losses[i]):
                print "Diff at", e_losses[i], " = (", log_cdfs_0[i], "-", process_cdfs_0[i], ") = ", diff
