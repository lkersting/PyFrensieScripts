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
#  Adjoint Electroatom Tests
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'H-Native' )
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
forward_file_name = datadir + data_list.get( 'electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
forward_data = Native.ElectronPhotonRelaxationDataContainer( forward_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()

###
###  Adjoint Electroatom/Electroatom Core Test Check
###
print "\n----- Electroatom Classes -----\n"
adjoint_brem_cs = adjoint_data.getAdjointBremsstrahlungElectronCrossSection()
adjoint_excitation_cs = adjoint_data.getAdjointAtomicExcitationCrossSection()
adjoint_analog_cs = adjoint_data.getAdjointTotalElasticCrossSection()
forward_inelastic_cs = adjoint_data.getForwardInelasticElectronCrossSection()

energies = [1e-5, 1e-3, 20.0]
for energy in energies:
  index = 0
  for i in range(0, adjoint_energy_grid.size ):
      if adjoint_energy_grid[i] <= energy:
          index = i

  energy_0 = adjoint_energy_grid[index]
  brem_cs_0 = adjoint_brem_cs[index]
  excitation_cs_0 = adjoint_excitation_cs[index]

  forward_inelastic_0 = forward_inelastic_cs[index]

  brem_cs = brem_cs_0
  excitation_cs = excitation_cs_0
  forward_inelastic = forward_inelastic_0
  tot_cs_0 = brem_cs + excitation_cs
  tot_cs = tot_cs_0

  if adjoint_energy_grid[index] != energy:
    energy_1 = adjoint_energy_grid[index+1]
    brem_cs_1 = adjoint_brem_cs[index+1]
    excitation_cs_1 = adjoint_excitation_cs[index+1]
    forward_inelastic_1 = forward_inelastic_cs[index+1]
    tot_cs_1 = brem_cs_1 + excitation_cs_1

    lin_interp = (energy - energy_0)/(energy_1 - energy_0)
    brem_cs = brem_cs_0 + (brem_cs_1 - brem_cs_0)*lin_interp
    excitation_cs = excitation_cs_0 + (excitation_cs_1 - excitation_cs_0)*lin_interp
    forward_inelastic = forward_inelastic_0 + (forward_inelastic_1 - forward_inelastic_0)*lin_interp
    tot_cs = tot_cs_0 + (tot_cs_1 - tot_cs_0)*lin_interp

  tot_cs = brem_cs + excitation_cs

  print "\nenergy = ", energy
  print '\tbrem_cs        = ','%.16e' % brem_cs
  print '\texcitation_cs  = ','%.16e' % excitation_cs
  print '\n\ttot_cs               = ','%.16e' % tot_cs
  print '\tforward_inelastic_cs = ','%.16e' % forward_inelastic
  print '\tindex = ', index

###
###  Adjoint Electroatom Native Factory/Electroatom Factory Test Check
###
print "\n----- Electroatom Factory Classes -----\n"
print "\n----- H -----\n"
tot_adjoint_elastic_cs = adjoint_data.getAdjointTotalElasticCrossSection()
reduced_cutoff_ratio = adjoint_data.getReducedCutoffCrossSectionRatios()
adjoint_cutoff_cs = adjoint_data.getAdjointCutoffElasticCrossSection()
adjoint_rutherford_cs = adjoint_data.getAdjointScreenedRutherfordElasticCrossSection()
adjoint_sr_index = adjoint_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
angular_energy_grid = adjoint_data.getAdjointElasticAngularEnergyGrid()
moment_cs_reduction = adjoint_data.getAdjointMomentPreservingCrossSectionReduction()
subshells = adjoint_data.getSubshells()
shell = subshells[0]
adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(shell)
adjoint_cutoff_dist = Collision.createLogLogLogCorrelatedCutoffElasticDistribution(adjoint_data, 0.9, 1e-7)
adjoint_mp_reaction = Collision.createLogLogLogCorrelatedMomentPreservingElasticReaction(adjoint_data, 0.9, 1e-7)
adjoint_hybrid_reaction = Collision.createLogLogLogCorrelatedHybridElasticReaction(adjoint_data, 0.9, 1e-7)
adjoint_hybrid_dist = Collision.createLogLogLogCorrelatedHybridElasticDistribution(adjoint_data, 0.9, 1e-7)
adjoint_cutoff_reaction = Collision.createLogLogLogCorrelatedCutoffElasticReaction(adjoint_data, 0.9, 1e-7)

energies = [1e-5, 1e-3, 20.0]
#energies = [1.0, 10.0, 20.0]
for energy in energies:
  index = 0
  for i in range(0, adjoint_energy_grid.size ):
      if adjoint_energy_grid[i] <= energy:
          index = i

  energy_0 = adjoint_energy_grid[index]
  brem_cs_0 = adjoint_brem_cs[index]
  excitation_cs_0 = adjoint_excitation_cs[index]
  ionization_cs_0 = adjoint_ionization_cs[index]
  tot_elastic_cs_0 = tot_adjoint_elastic_cs[index]
  cs_0 = brem_cs_0 + excitation_cs_0 + ionization_cs_0 + tot_elastic_cs_0

  cutoff_0 = adjoint_cutoff_dist.evaluateCutoffCrossSectionRatio( energy_0 )
  tot_cutoff_cs_0 = adjoint_cutoff_cs[index]
  cutoff_cs_0 = tot_cutoff_cs_0*cutoff_0

  mp_i = 0
  cross_section_reduction = 0
  for i in range(0, angular_energy_grid.size ):
    if angular_energy_grid[i] <= energy_0:
        mp_i = i
  if angular_energy_grid[mp_i] == energy_0:
    cross_section_reduction = moment_cs_reduction[mp_i]
  else:
    mp_log_interp = numpy.log(energy_0/angular_energy_grid[mp_i])/(angular_energy_grid[mp_i+1]/angular_energy_grid[mp_i])
    cross_section_reduction = moment_cs_reduction[mp_i]*pow((moment_cs_reduction[mp_i+1]/moment_cs_reduction[mp_i]),mp_log_interp )

  cutoff_cdf = adjoint_cutoff_dist.evaluateCutoffCrossSectionRatio( energy_0 )
  rutherford_cs_0 = tot_elastic_cs_0 - tot_cutoff_cs_0
  moment_cs_0 = cross_section_reduction*(rutherford_cs_0 + (1.0 - cutoff_cdf)*tot_cutoff_cs_0)

  brem_cs = brem_cs_0
  excitation_cs = excitation_cs_0
  ionization_cs = ionization_cs_0
  tot_elastic_cs = tot_elastic_cs_0
  tot_cs = cs_0

  tot_cutoff_cs = tot_cutoff_cs_0
  cutoff_cs = cutoff_cs_0
  moment_cs = moment_cs_0
  hybrid_cs_0 = moment_cs_0 + cutoff_cs_0
  hybrid_cs = hybrid_cs_0

  if adjoint_energy_grid[index] != energy:
    energy_1 = adjoint_energy_grid[index+1]
    brem_cs_1 = adjoint_brem_cs[index+1]
    excitation_cs_1 = adjoint_excitation_cs[index+1]
    ionization_cs_1 = adjoint_ionization_cs[index+1]
    tot_elastic_cs_1 = tot_adjoint_elastic_cs[index+1]
    cs_1 = brem_cs_1 + excitation_cs_1 + ionization_cs_1 + tot_elastic_cs_1

    cutoff_1 = adjoint_cutoff_dist.evaluateCutoffCrossSectionRatio( energy_1 )
    tot_cutoff_cs_1 = adjoint_cutoff_cs[index+1]
    cutoff_cs_1 = tot_cutoff_cs_1*cutoff_1

    mp_i = 0
    cross_section_reduction = 0
    for i in range(0, angular_energy_grid.size ):
      if angular_energy_grid[i] <= energy_1:
          mp_i = i
    if angular_energy_grid[mp_i] == energy_1:
      cross_section_reduction = moment_cs_reduction[mp_i]
    else:
      mp_log_interp = numpy.log(energy_1/angular_energy_grid[mp_i])/(angular_energy_grid[mp_i+1]/angular_energy_grid[mp_i])
      cross_section_reduction = moment_cs_reduction[mp_i]*pow((moment_cs_reduction[mp_i+1]/moment_cs_reduction[mp_i]),mp_log_interp )

    cutoff_cdf = adjoint_cutoff_dist.evaluateCutoffCrossSectionRatio( energy_1 )
    rutherford_cs_1 = tot_elastic_cs_1 - tot_cutoff_cs_1
    moment_cs_1 = cross_section_reduction*(rutherford_cs_1 + (1.0 - cutoff_cdf)*tot_cutoff_cs_1)

    lin_interp = (energy-energy_0)/(energy_1-energy_0)
    brem_cs = brem_cs_0 + (brem_cs_1 - brem_cs_0)*lin_interp
    excitation_cs = excitation_cs_0 + (excitation_cs_1 - excitation_cs_0)*lin_interp
    ionization_cs = ionization_cs_0 + (ionization_cs_1 - ionization_cs_0)*lin_interp
    tot_elastic_cs = tot_elastic_cs_0 + (tot_elastic_cs_1 - tot_elastic_cs_0)*lin_interp

    tot_cutoff_cs = tot_cutoff_cs_0 + (tot_cutoff_cs_1 - tot_cutoff_cs_0)*lin_interp

    cutoff_cs = cutoff_cs_0 + (cutoff_cs_1 - cutoff_cs_0)*lin_interp
    cutoff_cs = adjoint_cutoff_reaction.getCrossSection( energy )
    moment_cs = moment_cs_0 + (moment_cs_1 - moment_cs_0)*lin_interp
    moment_cs = adjoint_mp_reaction.getCrossSection( energy )

  hybrid_cs = cutoff_cs+moment_cs
  hybrid_cs = adjoint_hybrid_reaction.getCrossSection( energy )
  sr_cs = tot_elastic_cs - tot_cutoff_cs
  tot_cs = brem_cs + excitation_cs + ionization_cs + tot_elastic_cs
  tot_cs_with_cutoff = brem_cs + excitation_cs + ionization_cs + hybrid_cs

  max_excitation = (excitation_cs)/tot_cs
  max_brem = (excitation_cs+brem_cs)/tot_cs
  max_ionization = (ionization_cs+excitation_cs+brem_cs)/tot_cs
  max_elastic = (ionization_cs+tot_elastic_cs+excitation_cs+brem_cs)/tot_cs

  print "\nenergy = ", energy
  print '\tbrem_cs        = ','%.16e' % brem_cs
  print '\texcitation_cs  = ','%.16e' % excitation_cs
  print '\tionization_cs  = ','%.16e' % ionization_cs

  print '\n\tcutoff_cs      = ','%.16e' % tot_cutoff_cs
  print '\trutherford_cs  = ','%.16e' % sr_cs
  print '\tmoment_cs      = ','%.16e' % moment_cs
  print '\ttot_elastic_cs = ','%.16e' % tot_elastic_cs
  print '\ttot_cs         = ','%.16e' % tot_cs

  print '\n\thybrid_cs      = ','%.16e' % hybrid_cs
  print '\ttot_cs_with_cutoff = ','%.16e' % tot_cs_with_cutoff

  print '\n\tindex = ', index
  print '\treduced cutoff ratio at index','%.16e' % reduced_cutoff_ratio[index]

  print '\n\tmax excitation random number = ','%.16e' % max_excitation
  print '\tmax brem random number       = ','%.16e' % max_brem
  print '\tmax ionization random number = ','%.16e' % max_ionization
  print '\tmax elastic random number    = ','%.16e' % max_elastic

  sampling_ratio = tot_cutoff_cs/tot_elastic_cs
  print '\tcutoff sampling ratio = ','%.16e' % sampling_ratio



###
###  Adjoint Electroatom Native Factory/Electroatom Factory Test Check
###
print "\n----- Electroatom Factory Classes -----\n"
print "\n----- C -----\n"

data_list = cs_list.get( 'C-Native' )
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()
adjoint_brem_cs = adjoint_data.getAdjointBremsstrahlungElectronCrossSection()
adjoint_excitation_cs = adjoint_data.getAdjointAtomicExcitationCrossSection()
tot_adjoint_elastic_cs = adjoint_data.getAdjointTotalElasticCrossSection()
adjoint_cutoff_cs = adjoint_data.getAdjointCutoffElasticCrossSection()
reduced_cutoff_ratio = adjoint_data.getReducedCutoffCrossSectionRatios()
subshells = adjoint_data.getSubshells()
adjoint_rutherford_cs = adjoint_data.getAdjointScreenedRutherfordElasticCrossSection()
adjoint_sr_index = adjoint_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
adjoint_ionization_1_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[0])
adjoint_ionization_2_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[1])
adjoint_ionization_3_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[2])
adjoint_ionization_4_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[3])

energies = [1e-5, 1e-3, 20.0]
for energy in energies:
  index = 0
  for i in range(0, adjoint_energy_grid.size ):
      if adjoint_energy_grid[i] <= energy:
          index = i

  energy_0 = adjoint_energy_grid[index]
  brem_cs_0 = adjoint_brem_cs[index]
  excitation_cs_0 = adjoint_excitation_cs[index]
  ionization_1_cs_0 = adjoint_ionization_1_cs[index]
  ionization_2_cs_0 = adjoint_ionization_2_cs[index]
  ionization_3_cs_0 = adjoint_ionization_3_cs[index]
  ionization_4_cs_0 = adjoint_ionization_4_cs[index]
  tot_elastic_cs_0 = tot_adjoint_elastic_cs[index]
  cs_0 = brem_cs_0 + excitation_cs_0 + ionization_cs_0 + tot_elastic_cs_0

  cutoff_0 = reduced_cutoff_ratio[index]
  tot_cutoff_cs_0 = adjoint_cutoff_cs[index]
  cutoff_cs_0 = tot_cutoff_cs_0*cutoff_0
  moment_cs_0 = 0
  sr_cs_0 = 0.0
  if index >= adjoint_sr_index:
      sr_cs_0 = adjoint_rutherford_cs[index-adjoint_sr_index]

  brem_cs = brem_cs_0
  excitation_cs = excitation_cs_0
  ionization_1_cs = ionization_1_cs_0
  ionization_2_cs = ionization_2_cs_0
  ionization_3_cs = ionization_3_cs_0
  ionization_4_cs = ionization_4_cs_0
  tot_elastic_cs = tot_elastic_cs_0
  tot_cs = cs_0

  tot_cutoff_cs = tot_cutoff_cs_0
  sr_cs = sr_cs_0

  cutoff_cs = cutoff_cs_0
  moment_cs = moment_cs_0

  if adjoint_energy_grid[index] != energy:
    energy_1 = adjoint_energy_grid[index+1]
    brem_cs_1 = adjoint_brem_cs[index+1]
    excitation_cs_1 = adjoint_excitation_cs[index+1]
    ionization_1_cs_1 = adjoint_ionization_1_cs[index+1]
    ionization_2_cs_1 = adjoint_ionization_2_cs[index+1]
    ionization_3_cs_1 = adjoint_ionization_3_cs[index+1]
    ionization_4_cs_1 = adjoint_ionization_4_cs[index+1]
    tot_elastic_cs_1 = tot_adjoint_elastic_cs[index+1]
    cs_1 = brem_cs_1 + excitation_cs_1 + ionization_cs_1 + tot_elastic_cs_1

    cutoff_1 = reduced_cutoff_ratio[index+1]
    tot_cutoff_cs_1 = adjoint_cutoff_cs[index+1]
    cutoff_cs_1 = tot_cutoff_cs_1*cutoff_1
    moment_cs_1 = 0
    sr_cs_1 = 0.0
    if index+1 >= adjoint_sr_index:
        sr_cs_1 = adjoint_rutherford_cs[index+1-adjoint_sr_index]
    if sr_cs_0 > 0.0:
      sr_cs = sr_cs_0 + (sr_cs_1 - sr_cs_0)*lin_interp
    else:
      sr_cs = sr_cs_1

    lin_interp = (energy - energy_0)/(energy_1 - energy_0)
    brem_cs = brem_cs_0 + (brem_cs_1 - brem_cs_0)*lin_interp
    excitation_cs = excitation_cs_0 + (excitation_cs_1 - excitation_cs_0)*lin_interp
    ionization_1_cs = ionization_1_cs_0 + (ionization_1_cs_1 - ionization_1_cs_0)*lin_interp
    ionization_2_cs = ionization_2_cs_0 + (ionization_2_cs_1 - ionization_2_cs_0)*lin_interp
    ionization_3_cs = ionization_3_cs_0 + (ionization_3_cs_1 - ionization_3_cs_0)*lin_interp
    ionization_4_cs = ionization_4_cs_0 + (ionization_4_cs_1 - ionization_4_cs_0)*lin_interp
    tot_elastic_cs = tot_elastic_cs_0 + (tot_elastic_cs_1 - tot_elastic_cs_0)*lin_interp

    tot_cutoff_cs = tot_cutoff_cs_0 + (tot_cutoff_cs_1 - tot_cutoff_cs_0)*lin_interp

    cutoff_cs = cutoff_cs_0 + (cutoff_cs_1 - cutoff_cs_0)*lin_interp
    moment_cs = moment_cs_0 + (moment_cs_1 - moment_cs_0)*lin_interp

  ionization_cs = ionization_1_cs + ionization_2_cs + ionization_3_cs + ionization_4_cs
  hybrid_cs = cutoff_cs+moment_cs
  tot_cs = brem_cs + excitation_cs + ionization_cs + tot_elastic_cs
  tot_cs_with_cutoff = brem_cs + excitation_cs + ionization_cs + hybrid_cs

  max_excitation = (excitation_cs)/tot_cs
  max_brem = (excitation_cs+brem_cs)/tot_cs
  max_ionization = (ionization_cs+excitation_cs+brem_cs)/tot_cs
  max_elastic = (ionization_cs+tot_elastic_cs+excitation_cs+brem_cs)/tot_cs

  print "\nenergy = ", energy
  print '\tbrem_cs        = ','%.16e' % brem_cs
  print '\texcitation_cs  = ','%.16e' % excitation_cs
  print '\tionization_cs  = ','%.16e' % ionization_cs

  print '\n\tcutoff_cs     = ','%.16e' % tot_cutoff_cs
  print '\trutherford_cs  = ','%.16e' % sr_cs
  print '\tmoment_cs      = ','%.16e' % moment_cs
  print '\ttot_elastic_cs = ','%.16e' % tot_elastic_cs
  print '\ttot_cs         = ','%.16e' % tot_cs

  print '\n\thybrid_cs      = ','%.16e' % hybrid_cs
  print '\ttot_cs_with_cutoff = ','%.16e' % tot_cs_with_cutoff

  print '\n\tindex = ', index

  print '\n\tmax excitation random number = ','%.16e' % max_excitation
  print '\tmax brem random number       = ','%.16e' % max_brem
  print '\tmax ionization random number = ','%.16e' % max_ionization
  print '\tmax elastic random number    = ','%.16e' % max_elastic

  sampling_ratio = tot_cutoff_cs/tot_elastic_cs
  print '\tcutoff sampling ratio = ','%.16e' % sampling_ratio


###
###  Adjoint Electro Collision Handler Test Check
###
print "\n----- Electron Collision Handler -----\n"
print "\n----- Si -----\n"

data_list = cs_list.get( 'Si-Native' )
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()
adjoint_brem_cs = adjoint_data.getAdjointBremsstrahlungElectronCrossSection()
adjoint_excitation_cs = adjoint_data.getAdjointAtomicExcitationCrossSection()
tot_adjoint_elastic_cs = adjoint_data.getAdjointTotalElasticCrossSection()
adjoint_cutoff_cs = adjoint_data.getAdjointCutoffElasticCrossSection()
subshells = adjoint_data.getSubshells()
adjoint_ionization_1_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[0])
adjoint_ionization_2_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[1])
adjoint_ionization_3_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[2])
adjoint_ionization_4_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[3])
adjoint_ionization_5_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[4])
adjoint_ionization_6_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[5])
adjoint_ionization_7_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[6])
forward_inelastic_cs = adjoint_data.getForwardInelasticElectronCrossSection()

energies = [1e-3, 20.0, 1.55]
for energy in energies:
  index = 0
  for i in range(0, adjoint_energy_grid.size ):
      if adjoint_energy_grid[i] <= energy:
          index = i

  energy_0 = adjoint_energy_grid[index]
  brem_cs_0 = adjoint_brem_cs[index]
  excitation_cs_0 = adjoint_excitation_cs[index]
  ionization_1_cs_0 = adjoint_ionization_1_cs[index]
  ionization_2_cs_0 = adjoint_ionization_2_cs[index]
  ionization_3_cs_0 = adjoint_ionization_3_cs[index]
  ionization_4_cs_0 = adjoint_ionization_4_cs[index]
  ionization_5_cs_0 = adjoint_ionization_5_cs[index]
  ionization_6_cs_0 = adjoint_ionization_6_cs[index]
  ionization_7_cs_0 = adjoint_ionization_7_cs[index]
  tot_elastic_cs_0 = tot_adjoint_elastic_cs[index]
  cs_0 = brem_cs_0 + excitation_cs_0 + ionization_cs_0 + tot_elastic_cs_0

  forward_inelastic_0 = forward_inelastic_cs[index]
  tot_cutoff_cs_0 = adjoint_cutoff_cs[index]

  brem_cs = brem_cs_0
  excitation_cs = excitation_cs_0
  ionization_1_cs = ionization_1_cs_0
  ionization_2_cs = ionization_2_cs_0
  ionization_3_cs = ionization_3_cs_0
  ionization_4_cs = ionization_4_cs_0
  ionization_5_cs = ionization_5_cs_0
  ionization_6_cs = ionization_6_cs_0
  ionization_7_cs = ionization_7_cs_0
  tot_elastic_cs = tot_elastic_cs_0
  tot_cs = cs_0

  tot_cutoff_cs = tot_cutoff_cs_0
  forward_inelastic = forward_inelastic_0
  forward_cs_0 = forward_inelastic_0 + tot_elastic_cs_0
  forward_cs = forward_cs_0

  if adjoint_energy_grid[index] != energy:
    energy_1 = adjoint_energy_grid[index+1]
    brem_cs_1 = adjoint_brem_cs[index+1]
    excitation_cs_1 = adjoint_excitation_cs[index+1]
    ionization_1_cs_1 = adjoint_ionization_1_cs[index+1]
    ionization_2_cs_1 = adjoint_ionization_2_cs[index+1]
    ionization_3_cs_1 = adjoint_ionization_3_cs[index+1]
    ionization_4_cs_1 = adjoint_ionization_4_cs[index+1]
    ionization_5_cs_1 = adjoint_ionization_5_cs[index+1]
    ionization_6_cs_1 = adjoint_ionization_6_cs[index+1]
    ionization_7_cs_1 = adjoint_ionization_7_cs[index+1]
    tot_elastic_cs_1 = tot_adjoint_elastic_cs[index+1]
    cs_1 = brem_cs_1 + excitation_cs_1 + ionization_cs_1 + tot_elastic_cs_1

    tot_cutoff_cs_1 = adjoint_cutoff_cs[index+1]
    forward_inelastic_1 = forward_inelastic_cs[index+1]
    forward_cs_1 = forward_inelastic_1 + tot_elastic_cs_1

    lin_interp = (energy - energy_0)/(energy_1 - energy_0)
    brem_cs = brem_cs_0 + (brem_cs_1 - brem_cs_0)*lin_interp
    excitation_cs = excitation_cs_0 + (excitation_cs_1 - excitation_cs_0)*lin_interp
    ionization_1_cs = ionization_1_cs_0 + (ionization_1_cs_1 - ionization_1_cs_0)*lin_interp
    ionization_2_cs = ionization_2_cs_0 + (ionization_2_cs_1 - ionization_2_cs_0)*lin_interp
    ionization_3_cs = ionization_3_cs_0 + (ionization_3_cs_1 - ionization_3_cs_0)*lin_interp
    ionization_4_cs = ionization_4_cs_0 + (ionization_4_cs_1 - ionization_4_cs_0)*lin_interp
    ionization_5_cs = ionization_5_cs_0 + (ionization_5_cs_1 - ionization_5_cs_0)*lin_interp
    ionization_6_cs = ionization_6_cs_0 + (ionization_6_cs_1 - ionization_6_cs_0)*lin_interp
    ionization_7_cs = ionization_7_cs_0 + (ionization_7_cs_1 - ionization_7_cs_0)*lin_interp
    tot_elastic_cs = tot_elastic_cs_0 + (tot_elastic_cs_1 - tot_elastic_cs_0)*lin_interp

    tot_cutoff_cs = tot_cutoff_cs_0 + (tot_cutoff_cs_1 - tot_cutoff_cs_0)*lin_interp
    forward_inelastic = forward_inelastic_0 + (forward_inelastic_1 - forward_inelastic_0)*lin_interp
    forward_cs = forward_cs_0 + (forward_cs_1 - forward_cs_0)*lin_interp


  ionization_cs = ionization_1_cs + ionization_2_cs + ionization_3_cs + ionization_4_cs + ionization_5_cs + ionization_6_cs + ionization_7_cs
  tot_cs = brem_cs + excitation_cs + ionization_cs + tot_elastic_cs

  max_ionization = (ionization_cs)/tot_cs
  max_excitation = (ionization_cs+excitation_cs)/tot_cs
  max_brem = (ionization_cs+excitation_cs+brem_cs)/tot_cs
  max_elastic = (ionization_cs+tot_elastic_cs+excitation_cs+brem_cs)/tot_cs

  print "\nenergy = ", energy
  print '\tbrem_cs          = ','%.16e' % brem_cs
  print '\texcitation_cs    = ','%.16e' % excitation_cs
  print '\tionization_1_cs  = ','%.16e' % ionization_1_cs
  print '\tionization_7_cs  = ','%.16e' % ionization_7_cs

  print '\n\tcutoff_cs      = ','%.16e' % tot_cutoff_cs
  print '\ttot_elastic_cs = ','%.16e' % tot_elastic_cs
  print '\ttot_cs         = ','%.16e' % tot_cs
  print '\tforward_cs     = ','%.16e' % forward_cs

  print '\n\tindex = ', index

  print '\n\tmax ionization random number = ','%.16e' % max_ionization
  print '\tmax excitation random number = ','%.16e' % max_excitation
  print '\tmax brem random number       = ','%.16e' % max_brem
  print '\tmax elastic random number    = ','%.16e' % max_elastic

  weight = tot_cs/forward_cs
  print '\tweight_factor = ','%.16e' % weight