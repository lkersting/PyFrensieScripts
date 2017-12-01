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
data_list = cs_list.get( 'H-Native' )

# -------------------------------------------------------------------------- ##
#  Adjoint Elastic Data
# -------------------------------------------------------------------------- ##
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()


tot_adjoint_elastic_cs = adjoint_data.getAdjointTotalElasticCrossSection()
adjoint_cutoff_cs = adjoint_data.getAdjointCutoffElasticCrossSection()
reduced_cutoff_ratio = adjoint_data.getReducedCutoffCrossSectionRatios()
adjoint_screen_rutherford_cs = adjoint_data.getAdjointScreenedRutherfordElasticCrossSection()
adjoint_screen_rutherford_index = adjoint_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
mp_reduction = adjoint_data.getAdjointMomentPreservingCrossSectionReduction()

mp_reaction = Collision.createLogLogLogExactMomentPreservingElasticReaction( adjoint_data, 0.9, 1e-15)
hybrid_reaction = Collision.createLogLogLogExactHybridElasticReaction( adjoint_data, 0.9, 1e-15)
cutoff_dist = Collision.createLogLogLogExactCutoffElasticDistribution( adjoint_data, 0.999999, 1e-15)

###
###  Hybrid Distribution/Reaction Unit Test Check
###
energy = 1e-5
print "energy = ", energy
moment_cs = mp_reaction.getCrossSection( energy )
cutoff_cs = adjoint_cutoff_cs[0]
ratio = reduced_cutoff_ratio[0]
hrbrid_cs = cutoff_cs*ratio + moment_cs

print "discrete cross section = " '%.18e' % moment_cs
print "cutoff cross section = " '%.18e' % cutoff_cs
print "reduced cutoff ratio = " '%.18e' % ratio
print "hybrid cross section = " '%.18e' % hrbrid_cs

energy = 1e-3
print "\nenergy = ", energy
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] < energy:
        index = i

energy_0 = adjoint_energy_grid[index]
print "   energy_0 = ", energy_0
moment_cs_0 = mp_reaction.getCrossSection( energy_0 )
cutoff_cs_0 = adjoint_cutoff_cs[index]
ratio_0 = reduced_cutoff_ratio[index]
hybrid_cs_0 = cutoff_cs_0*ratio_0 + moment_cs_0
print "   hybrid_cs_0 = ", hybrid_cs_0, "\n"

energy_1 = adjoint_energy_grid[index+1]
print "   energy_1 = ", energy_1
moment_cs_1 = mp_reaction.getCrossSection( energy_1 )
cutoff_cs_1 = adjoint_cutoff_cs[index+1]
ratio_1 = reduced_cutoff_ratio[index+1]
hybrid_cs_1 = cutoff_cs_1*ratio_1 + moment_cs_1
print "   hybrid_cs_1 = ", hybrid_cs_1


cs = hybrid_cs_0 + ( hybrid_cs_1 - hybrid_cs_0 )*( energy - energy_0 )/( energy_1 - energy_0 )
cs = hybrid_reaction.getCrossSection( energy )
print "hybrid cross section = " '%.18e' % cs


energy = 20.0
print "\nenergy = ", energy
moment_cs = mp_reaction.getCrossSection( energy )
cutoff_cs = adjoint_cutoff_cs[adjoint_cutoff_cs.size-1]
ratio = reduced_cutoff_ratio[reduced_cutoff_ratio.size-1]
hrbrid_cs = cutoff_cs*ratio + moment_cs

print "discrete cross section = " '%.18e' % moment_cs
print "cutoff cross section = " '%.18e' % cutoff_cs
print "reduced cutoff ratio = " '%.18e' % ratio
print "hybrid cross section = " '%.18e' % hrbrid_cs


#adjoint_analog_dist = Collision.createHybridElasticDistribution(adjoint_data, 0.9)


#print '%.16e' % adjoint_analog_dist.evaluateCDF( 1.0e-3, 0.54 )
#print '%.16e' % adjoint_analog_dist.evaluateCDF( 1.0e-3, 0.9 )
#print '%.16e' % adjoint_analog_dist.evaluateCDF( 1.0e-3, 0.9239 )
#print '%.16e' % adjoint_analog_dist.evaluateCDF( 1.0e-3, 0.9788926224755288 )

#discrete_angles = h_adjoint_data.getAdjointMomentPreservingElasticDiscreteAngles(1e-3)
#discrete_weights = h_adjoint_data.getAdjointMomentPreservingElasticWeights(1e-3)
#discrete_energy_grid = h_adjoint_data.getAdjointElasticAngularEnergyGrid()
#print discrete_angles
#print discrete_weights

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

