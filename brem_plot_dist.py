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
#  Electroionization Data
# -------------------------------------------------------------------------- ##
# Possible Elements ['H-Native', 'Al-Native', 'Pb-Native']
elements = ['Pb-Native']
# Possible Interpolation Policies ["LogLogLog", "LinLinLin", "LinLinLog"]
interps = ["LogLogLog", "LinLinLin"]
# Possible Interpolation Schemes ["Correlated", "Exact"]
schemes = ["Correlated", "Exact"]
# Show Relative difference in schemes (True/False)
show_difference = True
# Possible energies [1e-2, 1e-1, 1.0, 15.7, 20.0]
energies = [15.7]
# Step length between plot points
step = 0.001
cdf_values = numpy.append(numpy.arange(0.0,1.0-step, step), 1.0-1e-15 )

random_numbers = ([])
for i in cdf_values:
  random_numbers = numpy.append(random_numbers,i)
  random_numbers = numpy.append(random_numbers,i)

plot_number = 1
for z in elements:
    print "\n----------------------------"
    print "-----", z, "Tests -----"
    print "----------------------------"
    data_list = cs_list.get( z )
    file_name = datadir + data_list.get( 'electroatomic_file_path' )
    native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )

    # Select distribution parameters
    for interp in interps:
      print "Interp = ", interp
      print "----------------------------"
      # Plot the given energy
      for energy in energies:

        title = 'Bremsstrahlung Energy Loss CDF at ' + str(energy)
        fig = plt.figure(num=plot_number, figsize=(10,5))
        if show_difference and len( schemes ) == 2:
          ax1 = plt.subplot2grid((2,1),(0, 0), colspan=5)
        else:
          plt.subplot2grid((1,1),(0, 0), colspan=5)
        plt.xlabel('Energy Loss')
        plt.ylabel('CDF')
        plt.title( title )
        # plt.xlim(0.85,1.0)
        # plt.ylim(0.02,0.05)

        # Plot all schemes on one graph
        samples = numpy.zeros(shape=(len( schemes ),len( cdf_values )))
        differences = numpy.zeros(shape=(2,len( cdf_values )))
        scheme_number = 0
        for scheme in schemes:

          # Create the distribution
          brem_dist = Collision.createLogLogLogCorrelatedBremsstrahlungDistribution(native_data, 1e-12)
          if interp == "LogLogLog":
            if scheme == "Exact":
              brem_dist = Collision.createLogLogLogExactBremsstrahlungDistribution(native_data, 1e-12)
            elif scheme == "Stochastic":
              brem_dist = Collision.createLogLogLogStochasticBremsstrahlungDistribution(native_data, 1e-12)
          elif interp == "LinLinLog":
            if scheme == "Correlated":
              brem_dist = Collision.createLinLinLogCorrelatedBremsstrahlungDistribution(native_data, 1e-12)
            elif scheme == "Exact":
              brem_dist = Collision.createLinLinLogExactBremsstrahlungDistribution(native_data, 1e-12)
            elif scheme == "Stochastic":
              brem_dist = Collision.createLinLinLogStochasticBremsstrahlungDistribution(native_data, 1e-12)
          elif interp == "LinLinLin":
            if scheme == "Correlated":
              brem_dist = Collision.createLinLinLinCorrelatedBremsstrahlungDistribution(native_data, 1e-12)
            elif scheme == "Exact":
              brem_dist = Collision.createLinLinLinExactBremsstrahlungDistribution(native_data, 1e-12)
            elif scheme == "Stochastic":
              brem_dist = Collision.createLinLinLinStochasticBremsstrahlungDistribution(native_data, 1e-12)

          Prng.RandomNumberGenerator.setFakeStream(random_numbers)
          angle = 0.0
          for i in range(0,len(cdf_values)):
            samples[scheme_number,i], angle = brem_dist.sample( energy )

          label = interp + " " + scheme
          plt.plot( samples[ scheme_number ], cdf_values, label=label)
          plt.xscale('log')
          plt.yscale('log')
          scheme_number = scheme_number + 1

        plt.legend( loc=4)

        if show_difference and len( schemes ) == 2:
          plt.subplot2grid((2,1),(1, 0), colspan=5, sharex=ax1)

          for i in range(0,len(cdf_values)):
            differences[0,i]= abs(samples[0,i] - samples[1,i])
            differences[1,i]= differences[0,i]/samples[0,i]

          plt.ylabel('Difference')
          plt.plot( samples[0], differences[0], label="Absolute Differences" )
          plt.plot( samples[0], differences[1], label="Rel. Absolute Differences" )
          plt.xscale('log')
          plt.legend( loc=2)

        plot_number = plot_number + 1

plt.show()