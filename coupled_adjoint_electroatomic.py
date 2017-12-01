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
adjoint_moment_cs = adjoint_data.getAdjointMomentPreservingCrossSection()
adjoint_moment_index = adjoint_data.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex()

adjoint_excitation_cs = adjoint_data.getAdjointAtomicExcitationCrossSection()
adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(1)
adjoint_brem_cs = adjoint_data.getAdjointBremsstrahlungElectronCrossSection()

###
###  Coupled Distribution/Reaction Unit Test Check
###

adjoint_dist = Collision.createCoupledElasticDistribution(adjoint_data, "LinLinLog", "Simplified Union", True, 1e-7)

energy = 6.625e1
cutoff_pdf = adjoint_dist.evaluatePDF( energy, 0.999999 )
print cutoff_pdf
max_cdf = adjoint_dist.evaluateCDF( energy, 0.1 )
print max_cdf
random_numbers = [0.1,0.5,1.0-1e-15]
Prng.RandomNumberGenerator.setFakeStream(random_numbers)
print "------------------------------------------------"
print "energy ", energy
print "------------------------------------------------"

E_out, mu = adjoint_dist.sample( energy )
print '1st angle cosine = ','%.16e' % mu

E_out, mu = adjoint_dist.sample( energy )
print '2nd angle cosine = ','%.16e' % mu

E_out, mu = adjoint_dist.sample( energy )
print '3rd angle cosine = ','%.16e' % mu


energy = 1e-4
max_cdf = adjoint_dist.evaluateCDF( energy, 1.0 )

random_numbers = [0.1,0.5,1.0-1e-15]
Prng.RandomNumberGenerator.setFakeStream(random_numbers)
print "------------------------------------------------"
print "energy ", energy
print "------------------------------------------------"

E_out, mu = adjoint_dist.sample( energy )
print '1st angle cosine = ','%.16e' % mu

E_out, mu = adjoint_dist.sample( energy )
print '2nd angle cosine = ','%.16e' % mu

E_out, mu = adjoint_dist.sample( energy )
print '3rd angle cosine = ','%.16e' % mu

