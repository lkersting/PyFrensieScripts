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
#  Electron Data
# -------------------------------------------------------------------------- ##
# Possible Elements ['H-Native', 'Al-Native', 'Pb-Native']
elements = ['Al-Native']
# Possible energies .005 .0093 .01 .011 .0134 .015 .0173 .02 .0252 .03 .04 .0415 .05 .06 .0621 .0818 and .102 MeV
energies = [.0173, .02]
# Possible Interpolation Policies ["LogLogLog", "LinLinLin", "LinLinLog"]
interps = ["LogLogLog"]
# Possible Reactions ["Elastic", "Bremsstrahlung", "Ionization", "Excitation"]
reactions = ["Bremsstrahlung"]
# Possible Interpolation Schemes ["Unit-base", "Unit-base Correlated", "Correlated"]
schemes = ["Unit-base", "Unit-base Correlated", "Correlated"]
# Show Relative difference in schemes (True/False)
show_difference = True
# # Step length between plot points
step = 0.01
length = int(1.0/step)
# Eval Tolerance
tol = 1e-12

plot_number = 0
for z in elements:
  print "\n----------------------------"
  print "-----", z, "Tests -----"
  print "----------------------------"
  data_list = cs_list.get(z)
  file_name = datadir + data_list.get('electroatomic_file_path')
  native_data = Native.ElectronPhotonRelaxationDataContainer(file_name)

  # Get energy grid
  energy_grid = native_data.getBremsstrahlungEnergyGrid()

  # Get the atomic shells (24 total)
  subshells = native_data.getSubshells()
  num_shells = len(subshells)

  # Get cross sections
  elastic_cross_sections = native_data.getTotalElasticCrossSection()
  brem_cross_sections = native_data.getBremsstrahlungCrossSection()
  brem_index = native_data.getBremsstrahlungCrossSectionThresholdEnergyIndex()
  excitation_cross_sections = native_data.getAtomicExcitationCrossSection()
  excitation_index = native_data.getAtomicExcitationCrossSectionThresholdEnergyIndex()

  ionization_cross_sections = []
  ionization_indexes = []
  for i in range(0, num_shells):
    shell = subshells[i]
    ionization_cross_sections.append(native_data.getElectroionizationCrossSection(shell))
    ionization_indexes.append(native_data.getElectroionizationCrossSectionThresholdEnergyIndex(shell))

  # Create distributions
  analog_dist = Collision.createLogLogLogUnitBaseCoupledElasticDistribution(native_data, "Two D Union", tol)
  excitation_dist = Collision.createAtomicExcitationDistribution(native_data)

  ionization_dists = []
  for i in range(0, num_shells):
    shell = subshells[i]
    binding_energy = native_data.getSubshellBindingEnergy(shell)
    ionization_dist_i = Collision.createLogLogLogUnitBaseElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
    ionization_dists.append(ionization_dist_i)


  # Plot the given energy
  for energy in energies:
    print "Energy = ", energy
    print "----------------------------"

    # Find energy in energy grid (lower bin index)
    index = 0
    for i in range(0, len(energy_grid)):
      if energy_grid[i] < energy:
        index = i


    for interp in interps:

      for reaction in reactions:

        if reaction == "Elastic":
          for n in range(0, len(schemes)):
            if interp == "LogLogLog":
              if schemes[n] == "Unit-base":
                brem_dist = Collision.createLogLogLogUnitBaseCoupledElasticDistribution(native_data, "Two D Union", tol)
              if schemes[n] == "Unit-base Correlated":
                brem_dist = Collision.createLogLogLogUnitBaseCorrelatedCoupledElasticDistribution(native_data, "Two D Union", tol)
              if schemes[n] == "Correlated":
                brem_dist = Collision.createLogLogLogCorrelatedCoupledElasticDistribution(native_data, "Two D Union", tol)
            elif interp == "LinLinLin":
              if schemes[n] == "Unit-base":
                brem_dist = Collision.createLinLinLinUnitBaseCoupledElasticDistribution(native_data, "Two D Union", tol)
              if schemes[n] == "Unit-base Correlated":
                brem_dist = Collision.createLinLinLinUnitBaseCorrelatedCoupledElasticDistribution(native_data, "Two D Union", tol)
              if schemes[n] == "Correlated":
                brem_dist = Collision.createLinLinLinCorrelatedCoupledElasticDistribution(native_data, "Two D Union", tol)


        if reaction == "Bremsstrahlung":

          # Choose the e_losses at which to evaluate the distribution
          e_losses = numpy.logspace(numpy.log10(1e-7), numpy.log10(energy-1e-8), num=length)

          pdfs = numpy.zeros(shape=(len(schemes), len(e_losses)))
          cdfs = numpy.zeros(shape=(len(schemes), len(e_losses)))
          samples = numpy.zeros(shape=(len(schemes), len(e_losses)))
          labels = []

          for n in range(0, len(schemes)):
            if interp == "LogLogLog":
              labels.append( schemes[n] + " - log" )
              if schemes[n] == "Unit-base":
                brem_dist = Collision.createLogLogLogUnitBaseBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Unit-base Correlated":
                brem_dist = Collision.createLogLogLogUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
                print brem_dist.evaluatePDF( energy, 1e-7)
              if schemes[n] == "Correlated":
                brem_dist = Collision.createLogLogLogCorrelatedBremsstrahlungDistribution(native_data, tol)
            elif interp == "LinLinLin":
              labels.append( schemes[n] + " - lin" )
              if schemes[n] == "Unit-base":
                brem_dist = Collision.createLinLinLinUnitBaseBremsstrahlungDistribution(native_data, tol)
              elif schemes[n] == "Unit-base Correlated":
                brem_dist = Collision.createLinLinLinUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
              elif schemes[n] == "Correlated":
                brem_dist = Collision.createLinLinLinCorrelatedBremsstrahlungDistribution(native_data, tol)

            for i in range(0, len(e_losses)):
              pdfs[n,i] = brem_dist.evaluatePDF(energy, e_losses[i])
              cdfs[n,i] = brem_dist.evaluateCDF(energy, e_losses[i])
              # print energy, e_losses[i], pdfs[n,i], cdfs[n,i]

          # Get lower and upper energy and the interpolation alpha
          energy_0 = energy_grid[index]
          energy_1 = energy_grid[index+1]
          lin_E_alpha = (energy - energy_0)/(energy_1 - energy_0)
          log_E_alpha = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

          # Get lower energy bin data
          e_losses_0 = native_data.getBremsstrahlungPhotonEnergy(energy_0)
          pdfs_0 = numpy.zeros(shape=len(e_losses_0))
          cdfs_0 = numpy.zeros(shape=len(e_losses_0))
          for i in range(0, len(e_losses_0)):
            pdfs_0[i] = brem_dist.evaluatePDF(energy_0, e_losses_0[i])
            cdfs_0[i] = brem_dist.evaluateCDF(energy_0, e_losses_0[i])

          # Get lower energy bin data
          e_losses_1 = native_data.getBremsstrahlungPhotonEnergy(energy_1)
          pdfs_1 = numpy.zeros(shape=len(e_losses_1))
          cdfs_1 = numpy.zeros(shape=len(e_losses_1))
          for i in range(0, len(e_losses_1)):
            pdfs_1[i] = brem_dist.evaluatePDF(energy_1, e_losses_1[i])
            cdfs_1[i] = brem_dist.evaluateCDF(energy_1, e_losses_1[i])


          title = 'Bremsstrahlung Energy Loss PDF at ' + str(energy) + ' MeV'
          plot_number = plot_number + 1

          fig = plt.figure(num=plot_number, figsize=(10, 5))
          if show_difference:
              plt.subplot2grid((2, 1), (0, 0), colspan=5)
          else:
              plt.subplot2grid((1, 1), (0, 0), colspan=5)
          plt.xlabel('Unit-base Energy Loss')
          plt.ylabel('PDF')
          plt.title(title)

          percent_values = (e_losses - e_losses[0])/(energy - e_losses[0])
          percent_values_0 = (e_losses_0 - e_losses_0[0])/(energy_0 - e_losses_0[0])
          percent_values_1 = (e_losses_1 - e_losses_1[0])/(energy_1 - e_losses_1[0])

          label0 = str(energy_0) +' MeV'
          label1 = str(energy_1) +' MeV'

          plt.plot(percent_values_0, pdfs_0, label=label0)
          plt.plot(percent_values_1, pdfs_1, label=label1)

          for n in range(0, len(schemes)):
            #print n, schemes[n], labels[n]
            #print pdfs[n]
            plt.plot(percent_values, pdfs[n], label=labels[n])

          plt.xscale('log')
          plt.yscale('log')
          plt.legend( loc=3)

          # Plot differences
          if len(schemes) > 1 and show_difference:
            # Plot differences in pdf
            plt.subplot2grid((2,1),(1, 0), colspan=5 )
            plt.ylabel('PDF Relative Difference')
            plt.xlabel('Unit-base Energy Loss')

            # Plot Differences
            rel_diff = numpy.zeros(shape=len(e_losses))
            for n in range(1, len(schemes)):
              diff_label = schemes[n] + ' vs ' + schemes[0]

              # Calculate difference between pdfs
              for i in range(0,length):
                rel_diff[i] = abs(pdfs[0,i] - pdfs[n,i])/pdfs[0,i]

              plt.plot( percent_values, rel_diff, label=diff_label)


            plt.xscale('log')
            plt.yscale('log')
            # plt.xlim(y_min,1.0)
            #plt.autoscale()
            #plt.xlim(cdfs[0],cdfs[len(cdfs)-1])
            #plt.ylim(y_min,y_max)
            plt.legend( loc=4)

          plot_number = plot_number + 1

          title = 'Bremsstrahlung Energy Loss CDF at ' + str(energy) + ' MeV'
          fig = plt.figure(num=plot_number, figsize=(10,5))
          if show_difference:
            plt.subplot2grid((2,1),(0, 0), colspan=5)
          else:
            plt.subplot2grid((1,1),(0, 0), colspan=5)
          plt.xlabel('Unit-base Energy Loss')
          plt.ylabel('CDF')
          plt.title( title )

          plt.plot( percent_values_0, cdfs_0, label=label0)
          plt.plot( percent_values_1, cdfs_1, label=label1)

          for n in range(0, len(schemes)):
            plt.plot(percent_values, cdfs[n], label=labels[n])

          plt.xscale('log')
          plt.yscale('log')
          # plt.xlim(x_min,1.0)
          plt.legend( loc=4)

          # Plot differences
          if len(schemes) > 1 and show_difference:
            # Plot differences in cdf
            plt.subplot2grid((2,1),(1, 0), colspan=5 )
            plt.ylabel('CDF Relative Difference')
            plt.xlabel('Unit-base Energy Loss')

            # Plot Differences
            rel_diff = numpy.zeros(shape=len(e_losses))
            for n in range(1, len(schemes)):
              diff_label = schemes[n] + ' vs ' + schemes[0]

              # Calculate difference between pdfs
              for i in range(0,length):
                rel_diff[i] = abs(cdfs[0,i] - cdfs[n,i])/cdfs[0,i]

              plt.plot( percent_values, rel_diff, label=diff_label)


            plt.xscale('log')
            plt.yscale('log')
            # plt.xlim(y_min,1.0)
            #plt.autoscale()
            #plt.xlim(cdfs[0],cdfs[len(cdfs)-1])
            #plt.ylim(y_min,y_max)
            plt.legend( loc=4)

        if reaction == "Ionization":
          for n in range(0, len(schemes)):
            if interp == "LogLogLog":
              if schemes[n] == "Unit-base":
                brem_dist = Collision.createLogLogLogUnitBaseBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Unit-base Correlated":
                brem_dist = Collision.createLogLogLogUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Correlated":
                brem_dist = Collision.createLogLogLogCorrelatedBremsstrahlungDistribution(native_data, tol)
            elif interp == "LinLinLin":
              if schemes[n] == "Unit-base":
                brem_dist = Collision.createLinLinLinUnitBaseBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Unit-base Correlated":
                brem_dist = Collision.createLinLinLinUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Correlated":
                brem_dist = Collision.createLinLinLinCorrelatedBremsstrahlungDistribution(native_data, tol)

        if reaction == "Excitation":
          for n in range(0, len(schemes)):
            if interp == "LogLogLog":
              if schemes[n] == "Unit-base":
                brem_dist = Collision.createLogLogLogUnitBaseBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Unit-base Correlated":
                brem_dist = Collision.createLogLogLogUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Correlated":
                brem_dist = Collision.createLogLogLogCorrelatedBremsstrahlungDistribution(native_data, tol)
            elif interp == "LinLinLin":
              if schemes[n] == "Unit-base":
                brem_dist = Collision.createLinLinLinUnitBaseBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Unit-base Correlated":
                brem_dist = Collision.createLinLinLinUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Correlated":
                brem_dist = Collision.createLinLinLinCorrelatedBremsstrahlungDistribution(native_data, tol)


plt.show()
