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

adjoint_excitation_cs = adjoint_data.getAdjointAtomicExcitationCrossSection()
adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(1)
adjoint_brem_cs = adjoint_data.getAdjointBremsstrahlungElectronCrossSection()

###
###  Cutoff Distribution/Reaction Unit Test Check
###
adjoint_cutoff_dist = Collision.createLogLogLogCorrelatedCutoffElasticDistribution(adjoint_data, 1.0, 1e-7)
cutoff_cdf = adjoint_cutoff_dist.evaluateCDF( 20.0, 0.9 )
ratio = reduced_cutoff_ratio[reduced_cutoff_ratio.size -1]
cutoff_cs = adjoint_cutoff_cs[reduced_cutoff_ratio.size -1]


print 'cutoff_cdf = ','%.16e' % cutoff_cdf
print 'ratio      = ','%.16e' % ratio
print 'cutoff_cs = ','%.16e' % cutoff_cs

cutoff_cdf = adjoint_cutoff_dist.evaluateCDF( 1e-3, 0.9 )
print 'cutoff_cdf = ','%.16e' % cutoff_cdf

cutoff_cdf = adjoint_cutoff_dist.evaluateCDF( 1e-5, 0.9 )
print 'cutoff_cdf = ','%.16e' % cutoff_cdf

# change cutoff
adjoint_cutoff_dist = Collision.createLogLogLogCorrelatedCutoffElasticDistribution(adjoint_data, 0.9, 1e-7)

cutoff_cdf = adjoint_cutoff_dist.evaluateCutoffCrossSectionRatio( 20.0 )
print 'cutoff_cdf = ','%.16e' % cutoff_cdf

cutoff_cdf = adjoint_cutoff_dist.evaluateCutoffCrossSectionRatio( 1e-3 )
print 'cutoff_cdf = ','%.16e' % cutoff_cdf

cutoff_cdf = adjoint_cutoff_dist.evaluateCutoffCrossSectionRatio( 1e-5 )
print 'cutoff_cdf = ','%.16e' % cutoff_cdf


