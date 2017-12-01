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

datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/frensie/src/packages/test_files/'
#datadir = '/home/lkersting/research/frensie-repos/lkersting/src/packages/test_files/'

source = Teuchos.FileInputSource( datadir + '/cross_sections.xml' )
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList( xml_obj )

h_data_list = cs_list.get( 'Al-Native' )
pb_data_list = cs_list.get( 'Pb-Native' )
h_adjoint_file_name = datadir + h_data_list.get( 'adjoint_electroatomic_file_path' )
h_native_file_name = datadir + h_data_list.get( 'electroatomic_file_path' )
pb_native_file_name = datadir + pb_data_list.get( 'electroatomic_file_path' )

h_native_data = Native.ElectronPhotonRelaxationDataContainer( h_native_file_name )
h_adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( h_adjoint_file_name )
pb_native_data = Native.ElectronPhotonRelaxationDataContainer( pb_native_file_name )

energy_grid = h_native_data.getElectronEnergyGrid()
angular_energy_grid = h_native_data.getElasticAngularEnergyGrid()
print angular_energy_grid
shells = h_native_data.getSubshells()

for shell in shells:
    print h_native_data.getSubshellBindingEnergy(shell)


angular_energy_grid = pb_native_data.getElasticAngularEnergyGrid()
print angular_energy_grid
shells = pb_native_data.getSubshells()

for shell in shells:
    print pb_native_data.getSubshellBindingEnergy(shell)

#print h_native_data.getMomentPreservingElasticDiscreteAngles( 1e-3 )
#lower = h_native_data.getMomentPreservingElasticDiscreteAngles( 8e-3 )
#upper = h_native_data.getMomentPreservingElasticDiscreteAngles( 1.6e-2 )

#print '%.16e'% lower[0]
#print '%.16e'% lower[1]
#print '%.16e'% upper[0]
#print '%.16e'% upper[1]
#print h_native_data.getMomentPreservingElasticDiscreteAngles( 1e5 )

#dist = Collision.createMomentPreservingElasticDistribution(h_native_data, 0.9, "LinLinLog", True, 1e-7 )

###
###  Forward Elastic Unit Test Check
###
#tot_elastic_cs = h_native_data.getTotalElasticCrossSection()
#cutoff_cs = h_native_data.getCutoffElasticCrossSection()
#screen_rutherford_cs = h_native_data.getScreenedRutherfordElasticCrossSection()
#screen_rutherford_index = h_native_data.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
#moment_cs = h_native_data.getMomentPreservingCrossSection()
#moment_index = h_native_data.getMomentPreservingCrossSectionThresholdEnergyIndex()

#dist     = Collision.createAnalogElasticDistribution(pb_native_data, True )
#lin_dist = Collision.createAnalogElasticDistribution(pb_native_data, False )

#print 'LOG: PDF', '%.16e'% dist.evaluatePDF(1e-4, -1.0)
#print 'LIN: PDF', '%.16e'% lin_dist.evaluatePDF(1e-4, -1.0)
## Sample
#energies = [1e-4, 1e-3,6.625e1]
#for i in energies:
#    print 'energy =',i
##    random_numbers = [2.23821564049245E-01, 2.23821564049245E-01,\
##                      3.60131793351356E-01, 3.60131793351356E-01,\
##                      0.4, 0.4,\
##                      0.453, 0.453]
#    random_numbers = [0.5, 0.5]
#    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#    for j in range(0,1):
#        e_out, pdf = dist.sample( i )
#        e_out, lin_pdf = lin_dist.sample( i )
#        print 'LOG: sample(','%.1e'% i,') = ','%.16e'% pdf,\
#            '\tLIN: sample(','%.3e'% i,') = ','%.16e'% lin_pdf

#dist = Collision.createHybridElasticDistribution(h_native_data, 0.9 )
#lin_dist = Collision.createHybridElasticDistribution(h_native_data, 0.9 )

#energies = [200.0]
#angles1 = [-0.01, 0.71, 0.999999, 1.0]
#angles2 = [0.9, 0.91, 9.23783127169921725e-01, 9.81773163837444063e-01]
#angles3 = [0.9, 9.33076191050687509e-01, 0.95, 9.99076800496474626e-01]
#energy = 1e-4
#angle = 0.9
#print lin_dist.sample( energy)
#for i in angles1:
#    energy = energies[0]
#    pdf = dist.evaluateCDF( energy, i )
#    lin_pdf = lin_dist.evaluateCDF( energy, i )
#    print 'LOG: PDF( E, mu ) = PDF(', energy,',','%.6e'% i,') = ','%.16e'% pdf,\
#          '\t LIN: PDF( E, mu ) = PDF(', energy,',','%.6e'% i,') = ','%.16e'% lin_pdf

#for i in angles2:
#    energy = energies[1]
#    pdf = dist.evaluate( energy, i )
#    lin_pdf = lin_dist.evaluateCDF( energy, i )
#    print 'LOG: PDF( E, mu ) = PDF(', energy,',','%.6e'% i,') = ','%.16e'% pdf,\
#          '\t LIN: PDF( E, mu ) = PDF(', energy,',','%.6e'% i,') = ','%.16e'% lin_pdf

#for i in angles3:
#    energy = energies[2]
#    pdf = dist.evaluate( energy, i )
#    lin_pdf = lin_dist.evaluateCDF( energy, i )
#    print 'LOG: PDF( E, mu ) = PDF(', energy,',','%.6e'% i,') = ','%.16e'% pdf,\
#          '\t LIN: PDF( E, mu ) = PDF(', energy,',','%.6e'% i,') = ','%.16e'% lin_pdf

#cutoff_cdf = dist.evaluateCDF( 1.99526E-04, 0.9 )
#print 'cutoff_cdf = ','%.16e' % cutoff_cdf
#cutoff_cdf = cut_dist.evaluateCDF( 1.99526E-04, 0.9 )
#print 'cutoff_cdf = ','%.16e' % cutoff_cdf

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == 1.99526E-04:
#        print 'energy = ', energy_grid[i]
#        print 'cutoff_cs = ', cutoff_cs[i]
#        print 'moment_preserving_cs = ','%.16e' % moment_cs[i-moment_index]
#        print 'hybrid_cs = ','%.16e' % (cutoff_cs[i]*cutoff_cdf + moment_cs[i-moment_index])

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == 4e-4:
#        print 'energy = ', energy_grid[i]
#        print 'moment_preserving_cs = ','%.16e' % moment_cs[i-moment_index]


#index_1 = 355
#index_2 = 365
#print moment_cs[index_1:index_2]
#print energy_grid[index_1:index_2]

#print 'moment preserving cs = ', '%.16e' % moment_cs[363]
#print 'energy = ', '%.16e' % energy_grid[363]

##
##  Brem Unit Test Check
##
#brem_lin = Collision.createBremsstrahlungDistribution(h_native_data, 82, True)
#print "Sample Brem"
#random_numbers = [0.5, 0.5, 0.49, 0.5, 0.48]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print brem_lin.sample( 1.0 )
#print '%.16e' % brem_lin.evaluate( 1.1e-5, 1.0 )
#print '%.16e' % brem_lin.evaluate( 20.0, 20.000000101 )
#print '%.16e' % brem_lin.evaluate( 21.0, 22.0 )
#print '%.16e' % brem_lin.evaluate( 1.0e-4, 1.0 )

###
###  Electroionization Unit Test Check
###
#subshells = pb_native_data.getSubshells()
#shell_0 = subshells[0]
#binding_energy = pb_native_data.getSubshellBindingEnergy( shell_0 )

##dist = Collision.createElectroionizationSubshellDistribution(pb_native_data, shell_0, binding_energy, True, False)
##lin_dist = Collision.createElectroionizationSubshellDistribution(pb_native_data, shell_0, binding_energy, False, False)

#LinLog = True
#LinLin = False
#Correlated = True
#Stochastic = False
#UnitBased = True
#Exact = False
#dist = Collision.createBremsstrahlungDistribution(pb_native_data, LinLog, Correlated, UnitBased)
#dist2 = Collision.createBremsstrahlungDistribution(pb_native_data, LinLog, Stochastic, Exact)

#brem_energy_grid = pb_native_data.getBremsstrahlungEnergyGrid()
#brem_dist1 = pb_native_data.getBremsstrahlungPhotonEnergy(2.0)
#brem_dist2 = pb_native_data.getBremsstrahlungPhotonPDF(2.0)
#print brem_energy_grid
#print brem_dist1
#print brem_dist2

#print dist.evaluatePDF( 2.1, 1e-7 )

#E = 1e-1
#E_1 = brem_energy_grid[4] # 1.18921e-02
#E_2 = brem_energy_grid[5] # 3.16228e-01

## LinLog
#E_p = E_1 + (E_2 - E_1)*numpy.log(E/E_1)/numpy.log(E_2/E_1)
#print E_p
## LinLin
#E_p = E_1 + (E_2 - E_1)*(E-E_1)/(E_2-E_1)
#print E_p

#num = 100
#random_numbers = [None]*(num+1)
#for i in range( 0, num):
#    random_numbers[i] = 0.01*i
#random_numbers[num] = 1.0-1e-10

#E = 1e2
#E_0 = brem_energy_grid[7]
#E_1 = brem_energy_grid[8]

#A = 0.0
#random_numbers = [A,A,A,A]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "------------------------------------------------"
#print "random number ", A
#print "------------------------------------------------"
#E_0_out, mu_0 = dist.sample( E_0 )
#E_1_out, mu_1 = dist.sample( E_1 )
#E_out, mu = dist.sample( E )
#print "E_0 = ", '%.6e' % E_0,": ", E_0_out
#print "E_1 = ", '%.6e' % E_1,": ", E_1_out
#print "E   = ", '%.6e' % E,": ", E_out
#print "------------------------------------------------"
#print "evaluateCDF(E,E_out) = ", A
#print "------------------------------------------------"
#print "E_0 = ", '%.6e' % E_0,": ", dist.evaluateCDF( E_0, E_0_out )
#print "E_1 = ", '%.6e' % E_1,": ", dist.evaluateCDF( E_1, E_1_out )
#print "E   = ", '%.6e' % E,": ", dist.evaluateCDF( E, E_out )
#print "------------------------------------------------\n"

#A = 0.5
#random_numbers = [A,A,A]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "------------------------------------------------"
#print "random number ", A
#print "------------------------------------------------"
#E_0_out, mu_0 = dist.sample( E_0 )
#E_1_out, mu_1 = dist.sample( E_1 )
#E_out, mu = dist.sample( E )
#print "E_0 = ", '%.6e' % E_0,": ", E_0_out
#print "E_1 = ", '%.6e' % E_1,": ", E_1_out
#print "E   = ", '%.6e' % E,": ", E_out
#print "------------------------------------------------"
#print "evaluateCDF(E,E_out) = ", A
#print "------------------------------------------------"
#print "E_0 = ", '%.6e' % E_0,": ", dist.evaluateCDF( E_0, E_0_out )
#print "E_1 = ", '%.6e' % E_1,": ", dist.evaluateCDF( E_1, E_1_out )
#print "E   = ", '%.6e' % E,": ", dist.evaluateCDF( E, E_out )
#print "------------------------------------------------\n"

#A = 1.0-1e-15
#random_numbers = [A,A,A]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "------------------------------------------------"
#print "random number ", A
#print "------------------------------------------------"
#E_0_out, mu_0 = dist.sample( E_0 )
#E_1_out, mu_1 = dist.sample( E_1 )
#E_out, mu = dist.sample( E )
#print "E_0 = ", '%.6e' % E_0,": ", E_0_out
#print "E_1 = ", '%.6e' % E_1,": ", E_1_out
#print "E   = ", '%.6e' % E,": ", E_out
#print "------------------------------------------------"
#print "evaluateCDF(E,E_out) = ", A
#print "------------------------------------------------"
#print "E_0 = ", '%.6e' % E_0,": ", dist.evaluateCDF( E_0, E_0_out )
#print "E_1 = ", '%.6e' % E_1,": ", dist.evaluateCDF( E_1, E_1_out )
#print "E   = ", '%.6e' % E,": ", dist.evaluateCDF( E, E_out )
#print "------------------------------------------------\n"

#A = dist.evaluateCDF( E, E_out )
#random_numbers = [A,A,A,A,A,A]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "------------------------------------------------"
#print "random number ", dist.evaluateCDF( E, E_out )
#print "------------------------------------------------"
#E_0_out, mu_0 = dist.sample( E_0 )
#E_1_out, mu_1 = dist.sample( E_1 )
#E_out, mu = dist.sample( E )
#print "E_0 = ", '%.6e' % E_0,": ", E_0_out
#print "E_1 = ", '%.6e' % E_1,": ", E_1_out
#print "E   = ", '%.6e' % E,": ", E_out

## LinLin
#cdf_0 = dist.evaluateCDF( E_0, E_0 )
#cdf_1 = dist.evaluateCDF( E_1, E_out )
#cdf = cdf_0 + (cdf_1 - cdf_0)*(E-E_0)/(E_1-E_0)
#print cdf


#for i in range( 0, num):
#    print dist.sample( 1e-1 )
#    angle, energy = dist.sample( 1e-1 )
#    print '%.16e' % energy
#print '%.16e' % dist.evaluatePDF( 9.12175e-2, 4.275e-4 )
#print '%.16e' % dist.evaluatePDF( 1e-1, 1e-2 )
#print '%.16e' % lin_dist.evaluatePDF( 1e-1, 1e-2 )



#h_adjoint_file_name = datadir + h_data_list.get( 'adjoint_electroatomic_file_path' )
#h_adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( h_adjoint_file_name )
#adjoint_energy_grid = h_adjoint_data.getAdjointElectronEnergyGrid()


##
##  Adjoint Elastic Unit Test Check
##
#tot_adjoint_elastic_cs = h_adjoint_data.getAdjointTotalElasticCrossSection()
#adjoint_cutoff_cs = h_adjoint_data.getAdjointCutoffElasticCrossSection()
#reduced_cutoff_ratio = h_adjoint_data.getReducedCutoffCrossSectionRatios()
#adjoint_screen_rutherford_cs = h_adjoint_data.getAdjointScreenedRutherfordElasticCrossSection()
#adjoint_screen_rutherford_index = h_adjoint_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
#adjoint_moment_cs = h_adjoint_data.getAdjointMomentPreservingCrossSection()
#adjoint_moment_index = h_adjoint_data.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex()

#adjoint_dist = Collision.createAnalogElasticDistribution(h_adjoint_data)
#energy = 1e-3
#energy_1 = adjoint_energy_grid[42]
#ratio_1 = reduced_cutoff_ratio[42]
#cutoff_1 = adjoint_cutoff_cs[42]
#moment_1 = adjoint_moment_cs[42]
#cross_section_1 = cutoff_1*ratio_1 + moment_1

#energy_2 = adjoint_energy_grid[43]
#ratio_2 = reduced_cutoff_ratio[43]
#cutoff_2 = adjoint_cutoff_cs[43]
#moment_2 = adjoint_moment_cs[43]
#cross_section_2 = cutoff_2*ratio_2 + moment_2

#linlin_slope = ( energy - energy_1 )/( energy_2 - energy_1 )
#ratio = ratio_1 + ( ratio_2 - ratio_1 )*linlin_slope
#cutoff = cutoff_1 + ( cutoff_2 - cutoff_1 )*linlin_slope
#moment = moment_1 + ( moment_2 - moment_1 )*linlin_slope
#cross_section = cross_section_1 + ( cross_section_2 - cross_section_1 )*linlin_slope

#print ' ratio = ', '%.16e' % ratio
#print ' cutoff = ', '%.16e' % cutoff
#print ' moment = ', '%.16e' % moment
#print ' cross section = ', '%.16e' % cross_section
 
