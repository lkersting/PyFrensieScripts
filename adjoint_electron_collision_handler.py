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
data_list = cs_list.get( 'Si-Native' )
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
forward_file_name = datadir + data_list.get( 'electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
forward_data = Native.ElectronPhotonRelaxationDataContainer( forward_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()

###
### getMacroscopicTotalCrossSection Test
###
adjoint_brem_cs = adjoint_data.getAdjointBremsstrahlungElectronCrossSection()
adjoint_excitation_cs = adjoint_data.getAdjointAtomicExcitationCrossSection()
adjoint_elastic_cs = adjoint_data.getAdjointTotalElasticCrossSection()
forward_inelastic_cs = adjoint_data.getForwardInelasticElectronCrossSection()

subshells = adjoint_data.getSubshells()
ionization_K_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[0])
ionization_L1_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[1])
ionization_L2_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[2])
ionization_L3_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[3])
ionization_M1_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[4])
ionization_M2_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[5])
ionization_M3_cs = adjoint_data.getAdjointElectroionizationCrossSection(subshells[6])

energies = [1e-3, 20.0, 1.55]
atomic_weight = 2.1441811932160583e-02
for energy in energies:
  # Get energy lower index
  index = 0
  for i in range(0, adjoint_energy_grid.size ):
      if adjoint_energy_grid[i] <= energy:
          index = i

  energy_0 = adjoint_energy_grid[index]
  brem_cs_0 = adjoint_brem_cs[index]
  excitation_cs_0 = adjoint_excitation_cs[index]
  elastic_cs_0 = adjoint_elastic_cs[index]
  K_cs_0 = ionization_K_cs[index]
  L1_cs_0 = ionization_L1_cs[index]
  L2_cs_0 = ionization_L2_cs[index]
  L3_cs_0 = ionization_L3_cs[index]
  M1_cs_0 = ionization_M1_cs[index]
  M2_cs_0 = ionization_M2_cs[index]
  M3_cs_0 = ionization_M3_cs[index]
  forward_inelastic_0 = forward_inelastic_cs[index]
  total_ionization_0 = K_cs_0 + L1_cs_0 + L2_cs_0 + L3_cs_0 + M1_cs_0 + M2_cs_0 + M3_cs_0
  cs_0 = brem_cs_0 + excitation_cs_0 + elastic_cs_0 + total_ionization_0

  if adjoint_energy_grid[index] != energy:
    energy_1 = adjoint_energy_grid[index+1]
    brem_cs_1 = adjoint_brem_cs[index+1]
    excitation_cs_1 = adjoint_excitation_cs[index+1]
    elastic_cs_1 = adjoint_elastic_cs[index+1]
    K_cs_1 = ionization_K_cs[index+1]
    L1_cs_1 = ionization_L1_cs[index+1]
    L2_cs_1 = ionization_L2_cs[index+1]
    L3_cs_1 = ionization_L3_cs[index+1]
    M1_cs_1 = ionization_M1_cs[index+1]
    M2_cs_1 = ionization_M2_cs[index+1]
    M3_cs_1 = ionization_M3_cs[index+1]
    forward_inelastic_1 = forward_inelastic_cs[index+1]
    total_ionization_1 = K_cs_1 + L1_cs_1 + L2_cs_1 + L3_cs_1 + M1_cs_1 + M2_cs_1 + M3_cs_1
    cs_1 = brem_cs_1 + excitation_cs_1 + elastic_cs_1 + total_ionization_1

    lin_slope = ( energy - energy_0 )/( energy_1 - energy_0 )

    brem_cs = (brem_cs_0 + (brem_cs_1 - brem_cs_0)*lin_slope)*atomic_weight
    excitation_cs = (excitation_cs_0 + (excitation_cs_1 - excitation_cs_0)*lin_slope)*atomic_weight
    elastic_cs = (elastic_cs_0 + (elastic_cs_1 - elastic_cs_0)*lin_slope)*atomic_weight
    K_cs = (K_cs_0 + (K_cs_1 - K_cs_0)*lin_slope)*atomic_weight
    L1_cs = (L1_cs_0 + (L1_cs_1 - L1_cs_0)*lin_slope)*atomic_weight
    L2_cs = (L2_cs_0 + (L2_cs_1 - L2_cs_0)*lin_slope)*atomic_weight
    L3_cs = (L3_cs_0 + (L3_cs_1 - L3_cs_0)*lin_slope)*atomic_weight
    M1_cs = (M1_cs_0 + (M1_cs_1 - M1_cs_0)*lin_slope)*atomic_weight
    M2_cs = (M2_cs_0 + (M2_cs_1 - M2_cs_0)*lin_slope)*atomic_weight
    M3_cs = (M3_cs_0 + (M3_cs_1 - M3_cs_0)*lin_slope)*atomic_weight

    forward_inelastic = (forward_inelastic_0 + (forward_inelastic_1 - forward_inelastic_0)*lin_slope)*atomic_weight
    total_ionization = (total_ionization_0 + (total_ionization_1 - total_ionization_0)*lin_slope)*atomic_weight
    cs = (cs_0 + (cs_1 - cs_0)*lin_slope)*atomic_weight
  else:
    brem_cs = brem_cs_0*atomic_weight
    excitation_cs = excitation_cs_0*atomic_weight
    elastic_cs = elastic_cs_0*atomic_weight
    K_cs = (K_cs_0)*atomic_weight
    L1_cs = (L1_cs_0)*atomic_weight
    L2_cs = (L2_cs_0)*atomic_weight
    L3_cs = (L3_cs_0)*atomic_weight
    M1_cs = (M1_cs_0)*atomic_weight
    M2_cs = (M2_cs_0)*atomic_weight
    M3_cs = (M3_cs_0)*atomic_weight
    forward_inelastic = forward_inelastic_0*atomic_weight
    cs = cs_0*atomic_weight
    total_ionization = total_ionization_0*atomic_weight
  forward_cs = forward_inelastic + elastic_cs
  weight_factor = cs/forward_cs

  print "\nenergy = ", energy
  print "\n----- Adjoint Cross Sections -----\n"
  print '\texcitation_cs    = ','%.16e' % excitation_cs
  print '\tbrem_cs          = ','%.16e' % brem_cs
  print '\telastic_cs       = ','%.16e' % elastic_cs
  print '\tionization_K_cs  = ','%.16e' % K_cs
  print '\tionization_M3_cs = ','%.16e' % M3_cs
  print '\t---------------------------------------'
  print '\ttot_adjoint_cs = ','%.16e' % cs
  print '\ttot_forward_cs = ','%.16e' % forward_cs

  ionization_ratio = total_ionization/cs
  brem_ratio = ( brem_cs )/cs
  excitation_ratio = ( excitation_cs )/cs
  elastic_ratio = ( elastic_cs )/cs

  print "\n----- Adjoint Cross Sections Ratios -----\n"
  print '\texcitation_cs ratio = ','%.16e' % excitation_ratio
  print '\tbrem_cs ratio       = ','%.16e' % brem_ratio
  print '\telastic_cs ratio    = ','%.16e' % elastic_ratio
  print '\tionization_cs ratio = ','%.16e' % ionization_ratio

  print "\n----- Adjoint Collide With Cell Material -----\n"

  max_ionization = total_ionization/cs
  max_brem = ( total_ionization + brem_cs )/cs
  max_excitation = ( total_ionization + brem_cs + excitation_cs )/cs
  max_elastic = ( total_ionization + brem_cs + excitation_cs + elastic_cs)/cs

  print '\tAdjoint_weight_factor = ','%.16e' % weight_factor
  print '\t---------------------------------------'
  print '\tmax ionization random number   = ','%.16e' % max_ionization
  print '\tmax brem random number         = ','%.16e' % max_brem
  print '\tactual min excitation random number   = .28545'
  print '\tactual max excitation random number   = ',excitation_cs/cs + .28545
  print '\tmax excitation random number   = ',max_excitation
  print '\tmax elastic random number      = ','%.16e' % max_elastic

  brem_dist = Collision.createLogLogLogCorrelatedBremsstrahlungDistribution(adjoint_data,1e-12)
  atomic_dist = Collision.createAtomicExcitationDistribution(adjoint_data)

  random_numbers = [ 0.0, 0.0, 0.0, 0.0,0.0, 0.0 ]
  Prng.RandomNumberGenerator.setFakeStream(random_numbers)

  incoming_energy = 1.55
  outgoing_energy, scattering_angle = brem_dist.sample( incoming_energy )
  print "\noutgoing_energy = ",'%.16e' % outgoing_energy
  print "scattering_angle = ",'%.16e' % scattering_angle

  outgoing_energy, scattering_angle = atomic_dist.sample( incoming_energy )
  print "\noutgoing_energy = ",'%.16e' % outgoing_energy
  print "scattering_angle = ",'%.16e' % scattering_angle
