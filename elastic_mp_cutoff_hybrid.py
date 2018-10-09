#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.DataGen.ElectronPhoton as EP
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.Utility.Interpolation as Interpolation
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
import numpy
import matplotlib.pyplot as plt

datadir = '/home/software/mcnpdata/'
#datadir = '/home/lkersting/frensie/src/packages/test_files/'

source = Teuchos.FileInputSource( datadir + '/cross_sections.xml' )
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList( xml_obj )

# -------------------------------------------------------------------------- ##
#  Native Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'H-Native' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
energy_grid = native_data.getElectronEnergyGrid()
elastic_energy_grid = native_data.getElasticAngularEnergyGrid()
cs_reduction = native_data.getMomentPreservingCrossSectionReduction()

hybrid_reaction = Collision.createHybridElasticReaction(native_data, 0.9, False, True, 1e-15)
cutoff_reaction = Collision.createCutoffElasticReaction(native_data, 1.0, False, True, 1e-15)
analog_reaction = Collision.createAnalogElasticReaction(native_data, False, True, 1e-15)

hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, False, True, 1e-15)
cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 1.0, False, True, 1e-15)
analog_dist = Collision.createAnalogElasticDistribution(native_data, False, True, 1e-15)

bin_index = 4
energy = elastic_energy_grid[bin_index]
mp_cs_reduction = cs_reduction[bin_index]

cutoff_cs = native_data.getCutoffElasticCrossSection()
moment_cs = native_data.getMomentPreservingCrossSection()
moment_index = native_data.getMomentPreservingCrossSectionThresholdEnergyIndex()

index = 0
for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

cut_cs = cutoff_cs[index]
mp_cs = moment_cs[index]
original_mp_cs = mp_cs/mp_cs_reduction

angles = native_data.getCutoffElasticAngles(energy)
pdfs = native_data.getCutoffElasticPDF(energy)
dcs = [None]*len(pdfs)
hybrid_dcs = [None]*len(pdfs)
hybrid_pdfs = [None]*len(pdfs)
hybrid_eval = [None]*len(pdfs)
hybrid_cs = hybrid_reaction.getCrossSection(energy)
cutoff_cdfs = [None]*len(pdfs)
for i in range(0, len(dcs)):
    dcs[i] = pdfs[i]*cut_cs
    hybrid_dcs[i] = hybrid_reaction.getDifferentialCrossSection(energy, angles[i] )
    hybrid_pdfs[i] = hybrid_dist.evaluatePDF(energy, angles[i] )
    hybrid_eval[i] = hybrid_dist.evaluate(energy, angles[i] )
    cutoff_cdfs[i] = cutoff_dist.evaluateCDF(energy, angles[i] )

size = 0
for i in range(0, angles.size ):
        if angles[i] <= 0.9:
            size = i+1

hybrid_angles = [None]*size
hybrid_cdf = [None]*size
for i in range(0, size-1):
    hybrid_angles[i] = angles[i]
    hybrid_cdf[i] = hybrid_dist.evaluateCDF(energy, angles[i] )

hybrid_angles[size-1] = 0.9
hybrid_cdf[size-1] = hybrid_dist.evaluateCDF(energy, hybrid_angles[size-1] )

discrete_angles = native_data.getMomentPreservingElasticDiscreteAngles(energy)
discrete_weights = native_data.getMomentPreservingElasticWeights(energy)

atomic_number = 1
min_electron_energy = 1e-5
max_electron_energy = 1e5
cutoff_angle_cosine = 0.9
tabular_evaluation_tol = 1e-15
linlinlog_interpolation_mode_on = False

hybrid_generator = EP.createMomentPreservingDataGenerator(
                atomic_number,
                native_data,
                min_electron_energy,
                max_electron_energy,
                cutoff_angle_cosine,
                tabular_evaluation_tol,
                linlinlog_interpolation_mode_on )

nodes, weights = hybrid_generator.evaluateDiscreteAnglesAndWeights( energy, 2 )

full_weights = [None] * 6
full_nodes = [None] * 6
full_nodes[0] = nodes[0]
full_nodes[1] = nodes[0]
full_nodes[2] = nodes[1]
full_nodes[3] = nodes[1]
full_nodes[4] = nodes[2]
full_nodes[5] = nodes[2]
full_weights[0] = 1e-12
full_weights[1] = weights[0]
full_weights[2] = weights[0]
full_weights[3] = weights[0] + weights[1]
full_weights[4] = weights[0] + weights[1]
full_weights[5] = weights[0] + weights[1] + weights[2]

cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, False, True, 1e-14)
cutoff_cdf = cutoff_dist.evaluateCutoffCrossSectionRatio( energy )

cross_section_ratio = cutoff_cs[index]*cutoff_cdf/moment_cs[index-moment_index]
sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
weight_0_ratio = (discrete_weights[0] + cross_section_ratio)/(1.0+cross_section_ratio)
weight_1_ratio = (discrete_weights[0] + discrete_weights[1] + cross_section_ratio)/(1.0+cross_section_ratio)


plot_weights = [None] * 5
plot_nodes = [None] * 5
plot_nodes[0] = cutoff_angle_cosine
plot_nodes[1] = discrete_angles[0]
plot_nodes[2] = discrete_angles[0]
plot_nodes[3] = discrete_angles[1]
plot_nodes[4] = discrete_angles[1]

plot_weights[0] = sampling_ratio
plot_weights[1] = sampling_ratio
plot_weights[2] = weight_0_ratio
plot_weights[3] = weight_0_ratio
plot_weights[4] = weight_1_ratio

fig1 = plt.figure(num=1, figsize=(10,5))
plt.xlabel('Angle Cosine', fontsize=15)
plt.ylabel('CDF', fontsize=15)
plt.title('Elastic Scattering Distributions', fontsize=16)
plt.xlim(0.7,1.005)
plt.ylim(5e-3,1.1)
plt.yscale('log')
plt.plot( angles, cutoff_cdfs, marker='o', linestyle='--', linewidth=2, label='Cutoff')
plt.plot( full_nodes, full_weights, marker='o', linestyle='--', linewidth=2, label='MP(0.9)')
plt.plot( hybrid_angles, hybrid_cdf, color='red', marker='o', linestyle='--', linewidth=2, label='Hybrid')
plt.plot( plot_nodes, plot_weights, color='red', marker='o', linestyle='--', linewidth=2)
plt.legend(loc=2)

fig1.savefig('./hybrid_cdf.pdf', bbox_inches='tight')

plt.show()

