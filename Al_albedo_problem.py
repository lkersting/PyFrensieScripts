#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
import numpy as numpy
import matplotlib.pyplot as plt

Utility.initFrensiePrng()

#datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/frensie/src/packages/test_files/'

source = Teuchos.FileInputSource(datadir + '/cross_sections.xml')
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList(xml_obj)

# -------------------------------------------------------------------------- ##
#  Electron Data
# -------------------------------------------------------------------------- ##
# Possible Elements ['H-Native', 'Al-Native', 'Pb-Native']
elements = ['Al-Native']
# Possible energies .005 .0093 .01 .011 .0134 .015 .0173 .02 .0252 .03 .04 .0415 .05 .06 .0621 .0818 and .102 MeV
energies = [.0173]
# Possible Interpolation Policies ["LogLogLog", "LinLinLin", "LinLinLog"]
interpolations = ["LogLogLog"]
# Possible Reactions ["Elastic", "Bremsstrahlung", "Ionization"]
reactions = ["Ionization"]
# Possible Interpolation Schemes ["Unit-base", "Unit-base Correlated", "Correlated"]
schemes = ["Unit-base", "Unit-base Correlated", "Correlated"]
# Possible comparisons ["PDF", "CDF", "Sample"]
comparisons = ["PDF"]
# Show Relative difference in schemes (True/False)
show_difference = True
# Just show first and last subshells
show_first_and_last_shell_only = True
# Step length between plot points
step = 0.001
length = int(1.0/step)
# Eval Tolerance
tol = 1e-10

plot_number = 0
dashes = ([8, 4, 2, 4], [8, 4], [6, 3, 3, 2, 3, 3], [4, 4])

for z in elements:
  print "\n----------------------------"
  print "-----", z, "Tests -----"
  print "----------------------------"
  data_list = cs_list.get(z)
  file_name = datadir + data_list.get('electroatomic_file_path')
  native_data = Native.ElectronPhotonRelaxationDataContainer(file_name)

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

  # Plot the given energy
  for energy in energies:
    print "Energy = ", energy
    print "----------------------------"

    for interpolation in interpolations:

      for reaction in reactions:

        if reaction == "Elastic":

          # Get energy grid
          energy_grid = native_data.getElasticAngularEnergyGrid()

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
          angles_0 = native_data.getCutoffElasticAngles(energy_0)

          # Get lower energy bin data
          angles_1 = native_data.getCutoffElasticAngles(energy_1)

          # Choose the angle cosines at which to evaluate the distribution
          second_angle = min(angles_0[1], angles_1[1])
          large_angles = numpy.linspace(second_angle, 0.1, num=length/10)
          large_angles = numpy.insert(large_angles, 0, -1.0)
          small_angles = numpy.logspace(numpy.log10(0.1), numpy.log10(1.0-1e-12), num=length*9/10)
          angles = numpy.unique(numpy.concatenate((large_angles, small_angles), axis=0))

          # Choose the random numbers at which to sample the distribution
          lower_number = numpy.logspace(numpy.log10(1e-10), numpy.log10(1e-4), num=length/10+1)
          upper_numbers = numpy.logspace(numpy.log10(1e-4), numpy.log10(1.0-1e-10), num=length*9/10)
          random_numbers = numpy.unique(numpy.concatenate((lower_number, upper_numbers), axis=0))

          pdfs = numpy.empty(shape=(len(schemes), length))
          cdfs = numpy.empty(shape=(len(schemes), length))
          samples = numpy.empty(shape=(len(schemes), length))

          labels = []

          for n in range(0, len(schemes)):
            if interpolation == "LogLogLog":
              labels.append(schemes[n] + " - log")
              if schemes[n] == "Unit-base":
                dist = Collision.createLogLogLogUnitBaseCoupledElasticDistribution(native_data, "One D Union", tol)
              if schemes[n] == "Unit-base Correlated":
                dist = Collision.createLogLogLogUnitBaseCorrelatedCoupledElasticDistribution(native_data, "Two D Union", tol)
              if schemes[n] == "Correlated":
                dist = Collision.createLogLogLogCorrelatedCoupledElasticDistribution(native_data, "Two D Union", tol)
            elif interpolation == "LinLinLin":
              labels.append(schemes[n] + " - lin")
              if schemes[n] == "Unit-base":
                dist = Collision.createLinLinLinUnitBaseCoupledElasticDistribution(native_data, "Two D Union", tol)
              if schemes[n] == "Unit-base Correlated":
                dist = Collision.createLinLinLinUnitBaseCorrelatedCoupledElasticDistribution(native_data, "Two D Union", tol)
              if schemes[n] == "Correlated":
                dist = Collision.createLinLinLinCorrelatedCoupledElasticDistribution(native_data, "Two D Union", tol)

            for i in range(0, length):
              pdfs[n, i] = dist.evaluatePDF(energy, angles[i])
              cdfs[n, i] = dist.evaluateCDF(energy, angles[i])

            if schemes[n] == "Unit-base":
              labels[n] = "Lower " + schemes[n] + " - log"
              upper_samples = numpy.zeros(len(random_numbers))
              Num = n

              # Sample both upper and lower distrbutions
              lower_random_numbers = numpy.zeros((len(random_numbers)*2))
              upper_random_numbers = numpy.zeros((len(random_numbers)*2))
              lower_random_numbers[0::2] = random_numbers
              upper_random_numbers[0::2] = random_numbers
              lower_random_numbers[1::2] = 1.0-1e-10

              # Sample lower distribution
              Prng.RandomNumberGenerator.setFakeStream(lower_random_numbers)
              for i in range(0, length):
                e_out, samples[n, i] = dist.sample(energy)

              # Sample upper distribution
              Prng.RandomNumberGenerator.setFakeStream(upper_random_numbers)
              for i in range(0, length):
                e_out, upper_samples[i] = dist.sample(energy)

            else:
              Prng.RandomNumberGenerator.setFakeStream(random_numbers)
              for i in range(0, length):
                e_out, samples[n, i] = dist.sample(energy)

          percent_values = (angles - angles[0])/(1.0 - angles[0])
          percent_values_0 = (angles_0 - angles_0[0])/(1.0 - angles_0[0])
          percent_values_1 = (angles_1 - angles_1[0])/(1.0 - angles_1[0])

          label0 = str(energy_0) +' MeV'
          label1 = str(energy_1) +' MeV'

          if "PDF" in comparisons:
            # Get lower energy bin data
            pdfs_0 = numpy.zeros(shape=len(angles_0))
            for i in range(0, len(angles_0)):
              pdfs_0[i] = dist.evaluatePDF(energy_0, angles_0[i])

            # Get lower energy bin data
            pdfs_1 = numpy.zeros(shape=len(angles_1))
            for i in range(0, len(angles_1)):
              pdfs_1[i] = dist.evaluatePDF(energy_1, angles_1[i])

            title = 'Elastic Angular PDF at ' + str(energy) + ' MeV'
            plot_number = plot_number + 1

            fig = plt.figure(num=plot_number, figsize=(10, 5))
            if len(schemes) > 1 and show_difference:
                plt.subplot2grid((2, 1), (0, 0), colspan=5)
            else:
                plt.subplot2grid((1, 1), (0, 0), colspan=5)
            plt.xlabel('Unit-base Angular Deflection')
            plt.ylabel('PDF')
            plt.title(title)

            plt.plot(percent_values_0, pdfs_0, label=label0)
            plt.plot(percent_values_1, pdfs_1, label=label1)

            for n in range(0, len(schemes)):
              plt.plot(percent_values, pdfs[n], dashes=dashes[n], label=labels[n])

            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc=2)

            # Plot differences
            if len(schemes) > 1 and show_difference:
              # Plot differences in pdf
              plt.subplot2grid((2, 1), (1, 0), colspan=5)
              plt.ylabel('PDF Relative Difference')
              plt.xlabel('Unit-base Angular Deflection')

              # Plot Differences
              rel_diff = numpy.zeros(shape=length)
              for n in range(1, len(schemes)):
                diff_label = schemes[0] + ' vs ' + schemes[n]

                # Calculate difference between pdfs
                for i in range(0, length):
                  rel_diff[i] = abs(pdfs[0, i] - pdfs[n, i])/pdfs[0, i]

                plt.plot(percent_values, rel_diff, dashes=dashes[n], label=diff_label)

              if len(schemes) > 2:
                for n in range(2, len(schemes)):
                  diff_label = schemes[1] + ' vs ' + schemes[n]

                  # Calculate difference between pdfs
                  rel_diff[0] = 0.0
                  rel_diff[length -1] = 0.0
                  for i in range(1, length-1):
                    rel_diff[i] = abs(pdfs[1, i] - pdfs[n, i])/pdfs[1, i]

                  plt.plot(percent_values, rel_diff, dashes=dashes[n], label=diff_label)

              plt.xscale('log')
              plt.yscale('log')
              plt.legend(loc=3)

          if "CDF" in comparisons:
            plot_number = plot_number + 1

            # Get lower energy bin data
            cdfs_0 = numpy.zeros(shape=len(angles_0))
            for i in range(0, len(angles_0)):
              cdfs_0[i] = dist.evaluateCDF(energy_0, angles_0[i])

            # Get lower energy bin data
            cdfs_1 = numpy.zeros(shape=len(angles_1))
            for i in range(0, len(angles_1)):
              cdfs_1[i] = dist.evaluateCDF(energy_1, angles_1[i])

            title = 'Elastic Angular CDF at ' + str(energy) + ' MeV'
            fig = plt.figure(num=plot_number, figsize=(10, 5))
            if len(schemes) > 1 and show_difference:
              plt.subplot2grid((2, 1), (0, 0), colspan=5)
            else:
              plt.subplot2grid((1, 1), (0, 0), colspan=5)
            plt.xlabel('Unit-base Angular Deflection')
            plt.ylabel('CDF')
            plt.title(title)

            plt.plot(percent_values_0, cdfs_0, label=label0)
            plt.plot(percent_values_1, cdfs_1, label=label1)

            for n in range(0, len(schemes)):
              plt.plot(percent_values, cdfs[n], dashes=dashes[n], label=labels[n])

            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc=2)

            # Plot differences
            if len(schemes) > 1 and show_difference:
              # Plot differences in cdf
              plt.subplot2grid((2, 1), (1, 0), colspan=5)
              plt.ylabel('CDF Relative Difference')
              plt.xlabel('Unit-base Angular Deflection')

              # Plot Differences
              rel_diff = numpy.zeros(shape=length)
              for n in range(1, len(schemes)):
                diff_label = schemes[0] + ' vs ' + schemes[n]

                # Calculate difference between cdfs
                rel_diff[0] = 0.0
                rel_diff[length -1] = 0.0
                for i in range(1, length-1):
                  rel_diff[i] = abs(cdfs[0, i] - cdfs[n, i])/cdfs[0, i]

                plt.plot(percent_values, rel_diff, dashes=dashes[n], label=diff_label)

              if len(schemes) > 2:
                for n in range(2, len(schemes)):
                  diff_label = schemes[1] + ' vs ' + schemes[n]

                  # Calculate difference between cdfs
                  rel_diff[0] = 0.0
                  rel_diff[length -1] = 0.0
                  for i in range(1, length-1):
                    rel_diff[i] = abs(cdfs[1, i] - cdfs[n, i])/cdfs[1, i]

                  plt.plot(percent_values, rel_diff, dashes=dashes[n-2], label=diff_label)

              plt.xscale('log')
              plt.yscale('log')
              plt.legend(loc=3)

          if "Sample" in comparisons:
            plot_number = plot_number + 1

            if "Unit-base" in schemes:
              samples = numpy.insert(samples, Num+1, upper_samples, axis=0)
              labels.insert(Num+1, "Upper Unit-base - log")
              schemes.insert(Num+1, "Upper Unit-base")
              schemes[Num] = "Lower Unit-base"

            # Get lower energy bin data
            samples_0 = numpy.zeros(shape=len(random_numbers))
            Prng.RandomNumberGenerator.setFakeStream(random_numbers)
            for i in range(0, len(random_numbers)):
              e_out, samples_0[i] = dist.sample(energy_0)

            # Get lower energy bin data
            samples_1 = numpy.zeros(shape=len(random_numbers))
            Prng.RandomNumberGenerator.setFakeStream(random_numbers)
            for i in range(0, len(random_numbers)):
              e_out, samples_1[i] = dist.sample(energy_1)

            percent_samples = numpy.zeros(shape=(len(schemes), length))
            for n in range(0, len(schemes)):
              percent_samples[n] = (samples[n] + 1.0)/(2.0)
            percent_samples_0 = (samples_0 + 1.0)/(2.0)
            percent_samples_1 = (samples_1 + 1.0)/(2.0)

            title = 'Sampled Elastic Angular Distribution at ' + str(energy) + ' MeV'
            fig = plt.figure(num=plot_number, figsize=(10, 5))
            if len(schemes) > 1 and show_difference:
              plt.subplot2grid((2, 1), (0, 0), colspan=5)
            else:
              plt.subplot2grid((1, 1), (0, 0), colspan=5)
            plt.ylabel('Unit-base Angular Deflection')
            plt.xlabel('CDF')
            plt.title(title)

            plt.plot(random_numbers, percent_samples_0, label=label0)
            plt.plot(random_numbers, percent_samples_1, label=label1)

            for n in range(0, len(schemes)):
              plt.plot(random_numbers, percent_samples[n], dashes=dashes[n], label=labels[n])

            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc=4)

            # Plot differences
            if len(schemes) > 1 and show_difference:
              # Plot differences in sampled values
              plt.subplot2grid((2, 1), (1, 0), colspan=5)
              plt.xlabel('CDF')
              plt.ylabel('Sampled Relative Difference')

              # Plot Differences
              rel_diff = numpy.zeros(shape=length)
              for m in range(0, len(schemes)-1):
                for n in range(m+1, len(schemes)):
                  diff_label = schemes[m] + ' vs ' + schemes[n]

                  # Calculate difference between sampled values
                  rel_diff[0] = 0.0
                  for i in range(0, length-1):
                    if samples[m, i] == 0.0:
                      rel_diff[i] = 0.0
                    else:
                      rel_diff[i] = abs((samples[m, i] - samples[n, i])/samples[m, i])

                  plt.plot(random_numbers, rel_diff, dashes=dashes[n], label=diff_label)

              plt.xscale('log')
              plt.yscale('log')
              plt.legend(loc=2)
              plt.autoscale()

        if reaction == "Bremsstrahlung":

          # Choose the e_losses at which to evaluate the distribution
          e_losses = numpy.logspace(numpy.log10(1e-7), numpy.log10(energy-1e-8), num=length)

          # Choose the random numbers at which to sample the distribution
          lower_number = numpy.logspace(numpy.log10(1e-10), numpy.log10(1e-4), num=length/10+1)
          upper_numbers = numpy.logspace(numpy.log10(1e-4), numpy.log10(1.0-1e-10), num=length*9/10)
          random_numbers = numpy.unique(numpy.concatenate((lower_number, upper_numbers), axis=0))
          angle_random_number = numpy.zeros((len(random_numbers)*2))
          angle_random_number[0::2] = random_numbers

          pdfs = numpy.zeros(shape=(len(schemes), length))
          cdfs = numpy.zeros(shape=(len(schemes), length))
          samples = numpy.zeros(shape=(len(schemes), length))
          labels = []

          Num = len(schemes)
          for n in range(0, len(schemes)):
            if interpolation == "LogLogLog":
              labels.append(schemes[n] + " - log")
              if schemes[n] == "Unit-base":
                dist = Collision.createLogLogLogUnitBaseBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Unit-base Correlated":
                dist = Collision.createLogLogLogUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
              if schemes[n] == "Correlated":
                dist = Collision.createLogLogLogCorrelatedBremsstrahlungDistribution(native_data, tol)
            elif interpolation == "LinLinLin":
              labels.append(schemes[n] + " - lin")
              if schemes[n] == "Unit-base":
                dist = Collision.createLinLinLinUnitBaseBremsstrahlungDistribution(native_data, tol)
              elif schemes[n] == "Unit-base Correlated":
                dist = Collision.createLinLinLinUnitBaseCorrelatedBremsstrahlungDistribution(native_data, tol)
              elif schemes[n] == "Correlated":
                dist = Collision.createLinLinLinCorrelatedBremsstrahlungDistribution(native_data, tol)

            if schemes[n] == "Unit-base":
              upper_samples = numpy.zeros(len(random_numbers))
              Num = n

              # Sample both upper and lower distrbutions
              lower_random_numbers = numpy.zeros((len(random_numbers)*3))
              upper_random_numbers = numpy.zeros((len(random_numbers)*3))
              lower_random_numbers[0::3] = random_numbers
              upper_random_numbers[0::3] = random_numbers
              lower_random_numbers[1::3] = 1.0-1e-10

              # Sample lower distribution
              Prng.RandomNumberGenerator.setFakeStream(lower_random_numbers)
              for i in range(0, length):
                samples[n, i], angle = dist.sample(energy)

              # Sample upper distribution
              Prng.RandomNumberGenerator.setFakeStream(upper_random_numbers)
              for i in range(0, length):
                upper_samples[i], angle = dist.sample(energy)

            else:
              Prng.RandomNumberGenerator.setFakeStream(angle_random_number)
              for i in range(0, length):
                samples[n, i], angle = dist.sample(energy)

            for i in range(0, length):
              pdfs[n, i] = dist.evaluatePDF(energy, e_losses[i])
              cdfs[n, i] = dist.evaluateCDF(energy, e_losses[i])

          # Get energy grid
          energy_grid = native_data.getBremsstrahlungEnergyGrid()

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

          e_losses_0 = native_data.getBremsstrahlungPhotonEnergy(energy_0)
          if e_losses[1] < e_losses_0[1]:
            e_losses_0 = numpy.insert(e_losses_0, 1, [e_losses[1]])

          e_losses_1 = native_data.getBremsstrahlungPhotonEnergy(energy_1)
          if e_losses[1] < e_losses_1[1]:
            e_losses_1 = numpy.insert(e_losses_1, 1, [e_losses[1]])

          if interpolation == "LinLinLin":
            percent_values = (e_losses - e_losses[0])/(energy - e_losses[0])
            percent_values_0 = (e_losses_0 - e_losses_0[0])/(energy_0 - e_losses_0[0])
            percent_values_1 = (e_losses_1 - e_losses_1[0])/(energy_1 - e_losses_1[0])
          elif interpolation == "LogLogLog":
            percent_values = numpy.log(e_losses/e_losses[0])/numpy.log(energy/e_losses[0])
            percent_values_0 = numpy.log(e_losses_0/e_losses_0[0])/numpy.log(energy_0/e_losses_0[0])
            percent_values_1 = numpy.log(e_losses_1/e_losses_1[0])/numpy.log(energy_1/e_losses_1[0])

          x_min = min(percent_values[1], min(percent_values_0[1], percent_values_1[1]))
          label0 = str(energy_0) +' MeV'
          label1 = str(energy_1) +' MeV'

          if "PDF" in comparisons:
            # Get lower energy bin data
            pdfs_0 = numpy.zeros(shape=len(e_losses_0))
            for i in range(0, len(e_losses_0)):
              pdfs_0[i] = dist.evaluatePDF(energy_0, e_losses_0[i])

            # Get lower energy bin data
            pdfs_1 = numpy.zeros(shape=len(e_losses_1))
            for i in range(0, len(e_losses_1)):
              pdfs_1[i] = dist.evaluatePDF(energy_1, e_losses_1[i])

            title = 'Bremsstrahlung Energy Loss PDF at ' + str(energy) + ' MeV'
            plot_number = plot_number + 1

            fig = plt.figure(num=plot_number, figsize=(10, 5))
            if len(schemes) > 1 and show_difference:
                plt.subplot2grid((2, 1), (0, 0), colspan=5)
            else:
                plt.subplot2grid((1, 1), (0, 0), colspan=5)
            plt.xlabel('Unit-base Energy Loss')
            plt.ylabel('PDF')
            plt.title(title)

            plt.plot(percent_values_0, pdfs_0, label=label0)
            plt.plot(percent_values_1, pdfs_1, label=label1)

            for n in range(0, len(schemes)):
              plt.plot(percent_values, pdfs[n], dashes=dashes[n], label=labels[n])

            plt.xscale('log')
            plt.yscale('log')
            plt.xlim(x_min, 1.0)
            plt.legend(loc=3)

            # Plot differences
            if len(schemes) > 1 and show_difference:
              # Plot differences in pdf
              plt.subplot2grid((2, 1), (1, 0), colspan=5)
              plt.ylabel('PDF Relative Difference')
              plt.xlabel('Unit-base Energy Loss')

              # Plot Differences
              rel_diff = numpy.zeros(shape=length)
              for m in range(0, len(schemes)-1):
                for n in range(m+1, len(schemes)):
                  diff_label = schemes[m] + ' vs ' + schemes[n]

                  # Calculate difference between pdfs
                  for i in range(0, length):
                    rel_diff[i] = abs(pdfs[m, i] - pdfs[n, i])/pdfs[m, i]

                  plt.plot(percent_values, rel_diff, dashes=dashes[n], label=diff_label)

              plt.xscale('log')
              plt.yscale('log')
              plt.xlim(x_min, 1.0)
              plt.legend(loc=4)

          if "CDF" in comparisons:
            plot_number = plot_number + 1

            # Get lower energy bin data
            cdfs_0 = numpy.zeros(shape=len(e_losses_0))
            for i in range(0, len(e_losses_0)):
              cdfs_0[i] = dist.evaluateCDF(energy_0, e_losses_0[i])

            # Get lower energy bin data
            cdfs_1 = numpy.zeros(shape=len(e_losses_1))
            for i in range(0, len(e_losses_1)):
              cdfs_1[i] = dist.evaluateCDF(energy_1, e_losses_1[i])

            title = 'Bremsstrahlung Energy Loss CDF at ' + str(energy) + ' MeV'
            fig = plt.figure(num=plot_number, figsize=(10, 5))
            if len(schemes) > 1 and show_difference:
              plt.subplot2grid((2, 1), (0, 0), colspan=5)
            else:
              plt.subplot2grid((1, 1), (0, 0), colspan=5)
            plt.xlabel('Unit-base Energy Loss')
            plt.ylabel('CDF')
            plt.title(title)

            plt.plot(percent_values_0, cdfs_0, label=label0)
            plt.plot(percent_values_1, cdfs_1, label=label1)

            for n in range(0, len(schemes)):
              plt.plot(percent_values, cdfs[n], dashes=dashes[n], label=labels[n])

            plt.xscale('log')
            plt.yscale('log')
            plt.xlim(x_min, 1.0)
            plt.legend(loc=4)

            # Plot differences
            if len(schemes) > 1 and show_difference:
              # Plot differences in cdf
              plt.subplot2grid((2, 1), (1, 0), colspan=5)
              plt.ylabel('CDF Relative Difference')
              plt.xlabel('Unit-base Energy Loss')

              # Plot Differences
              rel_diff = numpy.zeros(shape=length)
              for m in range(0, len(schemes)-1):
                for n in range(m+1, len(schemes)):
                  diff_label = schemes[m] + ' vs ' + schemes[n]

                  # Calculate difference between cdfs
                  rel_diff[0] = 0.0
                  rel_diff[length -1] = 0.0
                  for i in range(1, length-1):
                    rel_diff[i] = abs(cdfs[m, i] - cdfs[n, i])/cdfs[m, i]

                  plt.plot(percent_values, rel_diff, dashes=dashes[n], label=diff_label)

          if "Sample" in comparisons:
            plot_number = plot_number + 1

            if "Unit-base" in schemes:
              samples = numpy.insert(samples, Num+1, upper_samples, axis=0)
              labels[Num] = "Lower " + schemes[n] + " - log"
              labels.insert(Num+1, "Upper Unit-base - log")
              schemes.insert(Num+1, "Upper Unit-base")
              schemes[Num] = "Lower Unit-base"

            dist = Collision.createLogLogLogCorrelatedBremsstrahlungDistribution(native_data, tol)
            # Get lower energy bin data
            samples_0 = numpy.zeros(shape=len(random_numbers))
            Prng.RandomNumberGenerator.setFakeStream(angle_random_number)
            for i in range(0, len(random_numbers)):
              samples_0[i], angle = dist.sample(energy_0)

            # Get lower energy bin data
            samples_1 = numpy.zeros(shape=len(random_numbers))
            Prng.RandomNumberGenerator.setFakeStream(angle_random_number)
            for i in range(0, len(random_numbers)):
              samples_1[i], angle = dist.sample(energy_1)

            percent_samples = numpy.zeros(shape=(len(schemes), length))

            if interpolation == "LinLinLin":
              for n in range(0, len(schemes)):
                percent_samples[n] = (samples[n] - 1e-7)/(energy - 1e-7)
              percent_samples_0 = (samples_0 - 1e-7)/(energy_0 - 1e-7)
              percent_samples_1 = (samples_1 - 1e-7)/(energy_1 - 1e-7)
            elif interpolation == "LogLogLog":
              for n in range(0, len(schemes)):
                percent_samples[n] = numpy.log(samples[n]/1e-7)/numpy.log(energy/1e-7)
              percent_samples_0 = numpy.log(samples_0/1e-7)/numpy.log(energy_0/1e-7)
              percent_samples_1 = numpy.log(samples_1/1e-7)/numpy.log(energy_1/1e-7)

            title = 'Sampled Bremsstrahlung Energy Loss at ' + str(energy) + ' MeV'
            fig = plt.figure(num=plot_number, figsize=(10, 5))
            if len(schemes) > 1 and show_difference:
              plt.subplot2grid((2, 1), (0, 0), colspan=5)
            else:
              plt.subplot2grid((1, 1), (0, 0), colspan=5)
            plt.ylabel('Unit-base Energy Loss')
            plt.xlabel('CDF')
            plt.title(title)

            plt.plot(random_numbers, percent_samples_0, label=label0)
            plt.plot(random_numbers, percent_samples_1, label=label1)

            for n in range(0, len(schemes)):
              plt.plot(random_numbers, percent_samples[n], dashes=dashes[n], label=labels[n])

            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc=4)

            # Plot differences
            if len(schemes) > 1 and show_difference:
              # Plot differences in sampled values
              plt.subplot2grid((2, 1), (1, 0), colspan=5)
              plt.xlabel('CDF')
              plt.ylabel('Sampled Relative Difference')

              # Plot Differences
              rel_diff = numpy.zeros(shape=length)
              for m in range(0, len(schemes)-1):
                for n in range(m+1, len(schemes)):
                  diff_label = schemes[m] + ' vs ' + schemes[n]

                  # Calculate difference between sampled values
                  rel_diff[0] = 0.0
                  for i in range(0, length):
                    rel_diff[i] = abs((samples[m, i] - samples[n, i])/samples[m, i])

                  plt.plot(random_numbers, rel_diff, dashes=dashes[n], label=diff_label)

              plt.xscale('log')
              plt.yscale('log')
              plt.legend(loc=2)
              # plt.autoscale()

        if reaction == "Ionization":
          if show_first_and_last_shell_only:
            shells = [subshells[0], subshells[num_shells-1]]
          else:
            shells = subshells

          # Choose the random numbers at which to sample the distribution
          lower_number = numpy.logspace(numpy.log10(1e-10), numpy.log10(1e-4), num=length/10+1)
          upper_numbers = numpy.logspace(numpy.log10(1e-4), numpy.log10(1.0-1e-10), num=length*9/10)
          random_numbers = numpy.unique(numpy.concatenate((lower_number, upper_numbers), axis=0))

          schemes_original = list(schemes)
          for shell in shells:
            schemes = list(schemes_original)

            binding_energy = native_data.getSubshellBindingEnergy(shell)

            # Get energy grid
            energy_grid = native_data.getElectroionizationEnergyGrid(shell)

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
            e_losses_0 = native_data.getElectroionizationRecoilEnergy(shell, energy_0)

            # Get lower energy bin data
            e_losses_1 = native_data.getElectroionizationRecoilEnergy(shell, energy_1)

            # Choose the e_losses at which to evaluate the distribution
            max_energy = (energy - binding_energy)/2.0
            max_energy_interp = e_losses_0[len(e_losses_0)-1]*pow(e_losses_1[len(e_losses_1)-1]/e_losses_0[len(e_losses_0)-1], log_E_alpha)
            if interpolation == "LogLogLog":
              min_energy = e_losses_0[0]*pow(e_losses_1[0]/e_losses_0[0], log_E_alpha)
            elif interpolation == "LinLinLin":
              min_energy = energy_0 + (energy_1 - energy_0)*lin_E_alpha
            if min_energy > max_energy:
              print "No reaction is possible!"

            e_losses = numpy.logspace(numpy.log10(min_energy), numpy.log10(max_energy), num=length)
            if e_losses[1] < e_losses_0[1]:
              e_losses_0 = numpy.insert(e_losses_0, 1, [e_losses[1]])
            if e_losses[1] < e_losses_1[1]:
              e_losses_1 = numpy.insert(e_losses_1, 1, [e_losses[1]])

            pdfs = numpy.zeros(shape=(len(schemes), length))
            cdfs = numpy.zeros(shape=(len(schemes), length))
            samples = numpy.zeros(shape=(len(schemes), length))
            labels = []

            Num = len(schemes)
            for n in range(0, len(schemes)):
              if interpolation == "LogLogLog":
                labels.append(schemes[n] + " - log")
                if schemes[n] == "Unit-base":
                  dist = Collision.createLogLogLogUnitBaseElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
                if schemes[n] == "Unit-base Correlated":
                  dist = Collision.createLogLogLogUnitBaseCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
                if schemes[n] == "Correlated":
                  dist = Collision.createLogLogLogCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
              elif interpolation == "LinLinLin":
                labels.append(schemes[n] + " - lin")
                if schemes[n] == "Unit-base":
                  dist = Collision.createLinLinLinUnitBaseElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
                elif schemes[n] == "Unit-base Correlated":
                  dist = Collision.createLinLinLinUnitBaseCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
                elif schemes[n] == "Correlated":
                  dist = Collision.createLinLinLinCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)

              if schemes[n] == "Unit-base":
                upper_samples = numpy.zeros(len(random_numbers))
                Num = n

                # Sample both upper and lower distrbutions
                lower_random_numbers = numpy.zeros((len(random_numbers)*2))
                upper_random_numbers = numpy.zeros((len(random_numbers)*2))
                lower_random_numbers[0::2] = random_numbers
                upper_random_numbers[0::2] = random_numbers
                lower_random_numbers[1::2] = 1.0-1e-10

                # Sample lower distribution
                Prng.RandomNumberGenerator.setFakeStream(lower_random_numbers)
                for i in range(0, length):
                  samples[n, i], angle = dist.sample(energy)

                # Sample upper distribution
                Prng.RandomNumberGenerator.setFakeStream(upper_random_numbers)
                for i in range(0, length):
                  upper_samples[i], angle = dist.sample(energy)

              else:
                Prng.RandomNumberGenerator.setFakeStream(random_numbers)
                for i in range(0, length):
                  samples[n, i], angle = dist.sample(energy)

              for i in range(0, length-1):
                pdfs[n, i] = dist.evaluatePDF(energy, e_losses[i])
                cdfs[n, i] = dist.evaluateCDF(energy, e_losses[i])
              pdfs[n, length-1] = dist.evaluatePDF(energy, max_energy_interp)
              cdfs[n, length-1] = dist.evaluateCDF(energy, max_energy_interp)

            if interpolation == "LinLinLin":
              percent_values = (e_losses - e_losses[0])/(max_energy - e_losses[0])
              percent_values_0 = (e_losses_0 - e_losses_0[0])/(e_losses_0[len(e_losses_0)-1] - e_losses_0[0])
              percent_values_1 = (e_losses_1 - e_losses_1[0])/(e_losses_1[len(e_losses_1)-1] - e_losses_1[0])
            elif interpolation == "LogLogLog":
              percent_values = numpy.log(e_losses/e_losses[0])/numpy.log(max_energy/e_losses[0])
              percent_values_0 = numpy.log(e_losses_0/e_losses_0[0])/numpy.log(e_losses_0[len(e_losses_0)-1]/e_losses_0[0])
              percent_values_1 = numpy.log(e_losses_1/e_losses_1[0])/numpy.log(e_losses_1[len(e_losses_1)-1]/e_losses_1[0])

            x_min = min(percent_values[1], min(percent_values_0[1], percent_values_1[1]))
            label0 = str(energy_0) +' MeV'
            label1 = str(energy_1) +' MeV'

            if "PDF" in comparisons:
              # Get lower energy bin data
              pdfs_0 = numpy.zeros(shape=len(e_losses_0))
              for i in range(0, len(e_losses_0)):
                pdfs_0[i] = dist.evaluatePDF(energy_0, e_losses_0[i])

              # Get lower energy bin data
              pdfs_1 = numpy.zeros(shape=len(e_losses_1))
              for i in range(0, len(e_losses_1)):
                pdfs_1[i] = dist.evaluatePDF(energy_1, e_losses_1[i])

              title = 'Electro-ionization Energy Loss PDF at ' + str(energy) + ' MeV for shell ' + str(shell)
              plot_number = plot_number + 1

              fig = plt.figure(num=plot_number, figsize=(10, 5))
              if len(schemes) > 1 and show_difference:
                  plt.subplot2grid((2, 1), (0, 0), colspan=5)
              else:
                  plt.subplot2grid((1, 1), (0, 0), colspan=5)
              plt.xlabel('Unit-base Energy Loss')
              plt.ylabel('PDF')
              plt.title(title)

              plt.plot(percent_values_0, pdfs_0, label=label0)
              plt.plot(percent_values_1, pdfs_1, label=label1)

              for n in range(0, len(schemes)):
                plt.plot(percent_values, pdfs[n], dashes=dashes[n], label=labels[n])

              plt.xscale('log')
              plt.yscale('log')
              plt.xlim(x_min, 1.0)
              plt.legend(loc=3)

              # Plot differences
              if len(schemes) > 1 and show_difference:
                # Plot differences in pdf
                plt.subplot2grid((2, 1), (1, 0), colspan=5)
                plt.ylabel('PDF Relative Difference')
                plt.xlabel('Unit-base Energy Loss')

                # Plot Differences
                rel_diff = numpy.zeros(shape=length)
                for m in range(0, len(schemes)-1):
                  for n in range(m+1, len(schemes)):
                    diff_label = schemes[m] + ' vs ' + schemes[n]

                    # Calculate difference between pdfs
                    for i in range(0, length):
                      rel_diff[i] = abs(pdfs[m, i] - pdfs[n, i])/pdfs[m, i]

                    plt.plot(percent_values, rel_diff, dashes=dashes[n], label=diff_label)

                  plt.xscale('log')
                  plt.yscale('log')
                  plt.legend(loc=4)

            if "CDF"  in comparisons:
              plot_number = plot_number + 1

              # Get lower energy bin data
              cdfs_0 = numpy.zeros(shape=len(e_losses_0))
              for i in range(0, len(e_losses_0)):
                cdfs_0[i] = dist.evaluateCDF(energy_0, e_losses_0[i])

              # Get lower energy bin data
              cdfs_1 = numpy.zeros(shape=len(e_losses_1))
              for i in range(0, len(e_losses_1)):
                cdfs_1[i] = dist.evaluateCDF(energy_1, e_losses_1[i])

              title = 'Electro-ionization Energy Loss CDF at ' + str(energy) + ' MeV for shell ' + str(shell)
              fig = plt.figure(num=plot_number, figsize=(10, 5))
              if len(schemes) > 1 and show_difference:
                plt.subplot2grid((2, 1), (0, 0), colspan=5)
              else:
                plt.subplot2grid((1, 1), (0, 0), colspan=5)
              plt.xlabel('Unit-base Energy Loss')
              plt.ylabel('CDF')
              plt.title(title)

              plt.plot(percent_values_0, cdfs_0, label=label0)
              plt.plot(percent_values_1, cdfs_1, label=label1)

              for n in range(0, len(schemes)):
                plt.plot(percent_values, cdfs[n], dashes=dashes[n], label=labels[n])

              plt.xscale('log')
              plt.yscale('log')
              plt.xlim(x_min, 1.0)
              plt.legend(loc=4)

              # Plot differences
              if len(schemes) > 1 and show_difference:
                # Plot differences in cdf
                plt.subplot2grid((2, 1), (1, 0), colspan=5)
                plt.ylabel('CDF Relative Difference')
                plt.xlabel('Unit-base Energy Loss')

                # Plot Differences
                rel_diff = numpy.zeros(shape=length)
                for m in range(0, len(schemes)-1):
                  for n in range(m+1, len(schemes)):
                    diff_label = schemes[m] + ' vs ' + schemes[n]

                    # Calculate difference between cdfs
                    rel_diff[0] = 0.0
                    rel_diff[length -1] = 0.0
                    for i in range(1, length-1):
                      rel_diff[i] = abs(cdfs[m, i] - cdfs[n, i])/cdfs[m, i]

                    plt.plot(percent_values, rel_diff, dashes=dashes[n], label=diff_label)

                plt.xscale('log')
                plt.yscale('log')
                # plt.xlim(y_min,1.0)
                plt.legend(loc=4)

            if "Sample" in comparisons:
              plot_number = plot_number + 1

              if "Unit-base" in schemes:
                samples = numpy.insert(samples, Num+1, upper_samples, axis=0)
                labels[Num] = "Lower " + schemes[Num] + " - log"
                labels.insert(Num+1, "Upper Unit-base - log")
                schemes.insert(Num+1, "Upper Unit-base")
                schemes[Num] = "Lower Unit-base"

              dist = Collision.createLogLogLogCorrelatedElectroionizationSubshellDistribution(native_data, shell, binding_energy, tol)
              # Get lower energy bin data
              samples_0 = numpy.zeros(shape=len(random_numbers))
              Prng.RandomNumberGenerator.setFakeStream(random_numbers)
              for i in range(0, len(random_numbers)):
                samples_0[i], angle = dist.sample(energy_0)

              # Get lower energy bin data
              samples_1 = numpy.zeros(shape=len(random_numbers))
              Prng.RandomNumberGenerator.setFakeStream(random_numbers)
              for i in range(0, len(random_numbers)):
                samples_1[i], angle = dist.sample(energy_1)

              percent_samples = numpy.zeros(shape=(len(schemes), length))

              if interpolation == "LinLinLin":
                min_energy = e_losses_0[0] + (e_losses_1[0] - e_losses_0[0])*log_E_alpha

                for n in range(0, len(schemes)):
                  percent_samples[n] = (samples[n] - min_energy)/(max_energy - min_energy)
                percent_samples_0 = (samples_0 - e_losses_0[0])/(e_losses_0[len(e_losses_0)-1] - e_losses_0[0])
                percent_samples_1 = (samples_1 - e_losses_1[0])/(e_losses_1[len(e_losses_1)-1] - e_losses_1[0])
              elif interpolation == "LogLogLog":
                min_energy = e_losses_0[0]*pow(e_losses_1[0]/e_losses_0[0], log_E_alpha)

                for n in range(0, len(schemes)):
                  percent_samples[n] = numpy.log(samples[n]/min_energy)/numpy.log(max_energy/min_energy)
                percent_samples_0 = numpy.log(samples_0/e_losses_0[0])/numpy.log(e_losses_0[len(e_losses_0)-1]/e_losses_0[0])
                percent_samples_1 = numpy.log(samples_1/e_losses_1[0])/numpy.log(e_losses_1[len(e_losses_1)-1]/e_losses_1[0])


              title = 'Sampled Electro-ionization Energy Loss at ' + str(energy) + ' MeV for shell ' + str(shell)
              fig = plt.figure(num=plot_number, figsize=(10, 5))
              if len(schemes) > 1 and show_difference:
                plt.subplot2grid((2, 1), (0, 0), colspan=5)
              else:
                plt.subplot2grid((1, 1), (0, 0), colspan=5)
              plt.ylabel('Unit-base Angular Deflection')
              plt.xlabel('CDF')
              plt.title(title)

              plt.plot(random_numbers, percent_samples_0, label=label0)
              plt.plot(random_numbers, percent_samples_1, label=label1)

              for n in range(0, len(schemes)):
                plt.plot(random_numbers, percent_samples[n], dashes=dashes[n], label=labels[n])

              plt.xscale('log')
              plt.yscale('log')
              plt.legend(loc=4)

              # Plot differences
              if len(schemes) > 1 and show_difference:
                # Plot differences in sampled values
                plt.subplot2grid((2, 1), (1, 0), colspan=5)
                plt.xlabel('CDF')
                plt.ylabel('Sampled Relative Difference')

                # Plot Differences
                rel_diff = numpy.zeros(shape=length)
                for m in range(0, len(schemes)-1):
                  for n in range(m+1, len(schemes)):
                    diff_label = schemes[m] + ' vs ' + schemes[n]

                    # Calculate difference between sampled values
                    for i in range(0, length-1):
                      rel_diff[i] = abs((samples[m, i] - samples[n, i])/samples[m, i])

                    plt.plot(random_numbers, rel_diff, dashes=dashes[n], label=diff_label)

                plt.xscale('log')
                plt.yscale('log')
                plt.legend(loc=2)

plt.show()
