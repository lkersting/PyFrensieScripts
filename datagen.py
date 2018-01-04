#! /usr/bin/env python
import PyFrensie.Data.Native as Native
#import PyFrensie.DataGen.ElectronPhoton as EP
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
#import numpy
import matplotlib.pyplot as plt

Utility.initFrensiePrng()

### -------------------------------------------------------------------------- ##
###  Moment Preserving Reduction Check
### -------------------------------------------------------------------------- ##
# Get file names
file_1_angle = "./test_files/Au_MP_1_Angle.xml"
file_2_angle = "./test_files/Au_MP_2_Angles.xml"
file_4_angle = "./test_files/Au_MP_4_Angles.xml"
file_8_angle = "./test_files/Au_MP_8_Angles.xml"

# open files
mp_data_1 = Native.ElectronPhotonRelaxationDataContainer( file_1_angle )
mp_data_2 = Native.ElectronPhotonRelaxationDataContainer( file_2_angle )
mp_data_4 = Native.ElectronPhotonRelaxationDataContainer( file_4_angle )
mp_data_8 = Native.ElectronPhotonRelaxationDataContainer( file_8_angle )

energy_grid = mp_data_1.getElasticAngularEnergyGrid()
reductions_1 = mp_data_1.getMomentPreservingCrossSectionReduction()
reductions_2 = mp_data_2.getMomentPreservingCrossSectionReduction()
reductions_4 = mp_data_4.getMomentPreservingCrossSectionReduction()
reductions_8 = mp_data_8.getMomentPreservingCrossSectionReduction()

for i in range(0,len(reductions_1)):
    reductions_1[i] = 1.0/reductions_1[i]
    reductions_2[i] = 1.0/reductions_2[i]
    reductions_4[i] = 1.0/reductions_4[i]
    reductions_8[i] = 1.0/reductions_8[i]

fig1 = plt.figure(num=1, figsize=(10,5))
plt.subplot2grid((1,6),(0, 0), colspan=5)
plt.xlabel('E (MeV)')
plt.ylabel('Ratio of Total Cross Sections')
plt.title('Ratio of Total Cross Sections (Analog/MP) for H')
plt.xlim(1e-5,1e2)
plt.ylim(1.0,1e5)
plt.loglog( energy_grid, reductions_1, marker='o', linestyle='--', label='1-Angle')
plt.loglog( energy_grid, reductions_2, marker='o', linestyle='--', label='2-Angles')
plt.plot( energy_grid, reductions_4, marker='o', linestyle='--', label='4-Angles')
plt.plot( energy_grid, reductions_8, marker='o', linestyle='--', label='8-Angles')
plt.legend(loc=2)

fig1.savefig('./discrete_cs_ratios.pdf', bbox_inches='tight')

plt.show()

#atomic_number = 79
#min_electron_energy = 1e-5
#max_electron_energy = 1e5
#cutoff_angle_cosine = 0.9
#tabular_evaluation_tol = 1e-12
#linlinlog_interpolation_mode_on = True

#hybrid_generator = EP.createMomentPreservingDataGenerator(
#                atomic_number,
#                native_data,
#                min_electron_energy,
#                max_electron_energy,
#                cutoff_angle_cosine,
#                tabular_evaluation_tol,
#                linlinlog_interpolation_mode_on )

#cutoff_angle_cosine = -1.0

#mp_generator = EP.createMomentPreservingDataGenerator(
#                atomic_number,
#                native_data,
#                min_electron_energy,
#                max_electron_energy,
#                cutoff_angle_cosine,
#                tabular_evaluation_tol,
#                linlinlog_interpolation_mode_on )

#nodes, weights = hybrid_generator.evaluateDiscreteAnglesAndWeights( 1e5, 2 )
#print weights

#angles = 2
#for e in energy_grid:
#    print "energy =",e
#    nodes, weights = mp_generator.evaluateDiscreteAnglesAndWeights( e, angles )
#    print "\tweights =",weights
#    ratio = 0.0
#    for i in range(0, angles):
#        ratio += weights[i]
#    print ratio
#    ratio = 1.0 - weights[len(weights)-1]
#    print ratio



