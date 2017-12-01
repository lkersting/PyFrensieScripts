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
# Electroatom Tests
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'Pb-Native' )
native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )
energy_grid = native_data.getElectronEnergyGrid()
subshells = native_data.getSubshells()

###
### Electroatom/Electroatom Core Test Check
###
print "\n----- Electroatom Native Factory Class -----"
brem_cross_sections = native_data.getBremsstrahlungCrossSection()
brem_index = native_data.getBremsstrahlungCrossSectionThresholdEnergyIndex()
excitation_cross_sections = native_data.getAtomicExcitationCrossSection()

tot_elastic_cross_sections = native_data.getTotalElasticCrossSection()
cutoff_cross_sections = native_data.getCutoffElasticCrossSection()
moment_cross_sections = Collision.createLogLogLogExactMomentPreservingElasticReaction(native_data, 0.9, 1e-15)
hybrid_dist = Collision.createLogLogLogExactHybridElasticReaction(native_data, 0.9, 1e-15)
analog_dist = Collision.createLogLogLogExactCoupledElasticDistribution(native_data, "Two D Union",1e-15)
cutoff_dist = Collision.createLogLogLogExactCutoffElasticDistribution(native_data, 0.9, 1e-15)

energies = [1e-5, 2e-1, 1e5, 1e-3, 4e-4]
for energy in energies:

    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    print "\nEnergy = ",energy,'\tindex = ', index

    brem_cs = brem_cross_sections[index-brem_index]
    excitation_cs = excitation_cross_sections[index]
    cutoff_cs = cutoff_cross_sections[index]
    cutoff_cdf = analog_dist.evaluateCDFAtCutoff( energy )
    analog_cs = cutoff_cs/cutoff_cdf
    cutoff_ratio = cutoff_dist.evaluateCutoffCrossSectionRatio( energy )
    moment_cs = moment_cross_sections.getCrossSection( energy )
    hybrid_cs = cutoff_cs*cutoff_ratio + moment_cs
    hybrid_cs =  hybrid_dist.getCrossSection( energy )
    ionization_cs = 0.0

    for shell in subshells:
        ionization_cross_section = native_data.getElectroionizationCrossSection(shell)
        ionization_index = native_data.getElectroionizationCrossSectionThresholdEnergyIndex(shell)

        i = index-ionization_index
        if i >= 0:
            ionization_cs += ionization_cross_section[i]

    inelastic_cs = brem_cs + excitation_cs + ionization_cs
    total_cs_analog = inelastic_cs + analog_cs
    total_cs_cutoff = inelastic_cs + hybrid_cs

    print '\tbrem_cs       = ','%.16e' % brem_cs
    print '\texcitation_cs = ','%.16e' % excitation_cs
    print '\tionization_cs = ','%.16e' % ionization_cs
    print '\t------------------------------------------------'
    print '\tinelastic_cs  = ','%.16e' % inelastic_cs
    print '\t------------------------------------------------'
    print '\tanalog_cs (lin)  = ','%.16e' % analog_cs
    print '\tcutoff_cs (log) = ','%.16e' % cutoff_cs
    print '\tmoment_cs (log) = ','%.16e' % moment_cs
    print '\thybrid_cs (log) = ','%.16e' % hybrid_cs
    print '\t------------------------------------------------'
    print '\ttotal cs (analog) = ','%.16e' % total_cs_analog
    print '\ttotal cs (hybrid) = ','%.16e' % total_cs_cutoff

energies = [1.0e5, 1.995260e1, 6.309570e0, 1.995260e-3, 1.995260e-4, 1.0e-5]
print "\n--- Analog Cross Section ---"
for energy in energies:
    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    energy_0 = energy_grid[index]
    elastic_cs = 0.0
    if energy_0 != energy:
        energy_1 = energy_grid[index+1]
        lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )

        cutoff_cdf = analog_dist.evaluateCDFAtCutoff( energy_0 )
        elastic_cs_0 = cutoff_cross_sections[index]/cutoff_cdf

        cutoff_cdf = analog_dist.evaluateCDFAtCutoff( energy_1 )
        elastic_cs_1 = cutoff_cross_sections[index+1]/cutoff_cdf

        elastic_cs = elastic_cs_0 + (elastic_cs_1 - elastic_cs_0)*lin_interp
    else:
        cutoff_cdf = analog_dist.evaluateCDFAtCutoff( energy )
        elastic_cs = cutoff_cross_sections[index]/cutoff_cdf
    print '\tcs[','%.6e' %energy,']:','%.16e' % elastic_cs

energies = [1.0e5, 1.995260e1, 6.30957, 1e-3, 1.995260e-4, 1.0e-5]
print "\n--- Hybrid Cross Section ---"
analog_dist = Collision.createLogLogLogExactCutoffElasticDistribution(native_data, 1.0, 1e-15)
for energy in energies:
    print "\nEnergy = ",energy
    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    cutoff_cs = cutoff_cross_sections[index]
    cutoff_ratio = cutoff_dist.evaluateCutoffCrossSectionRatio( energy )
    print '\n','%.16e' % cutoff_ratio
    print '%.16e' % analog_dist.evaluateCDF( energy, 0.9 )
    moment_cs = moment_cross_sections.getCrossSection( energy )
    hybrid_cs = cutoff_cs*cutoff_ratio + moment_cs

    print '\tcutoff_cs','%.16e' % cutoff_cs
    print '\tmoment_cs','%.16e' % moment_cs
    print '\thybrid_cs','%.16e' % hybrid_cs


print "\n--- Electroatom Factory ---"
energies = [2e-3, 4e-4, 9e-5]
for energy in energies:

   index = 0
   for i in range(0, energy_grid.size ):
       if energy_grid[i] == energy:
           index = i

   print "\nEnergy = ",energy,'\tindex = ', index

   brem_cs = brem_cross_sections[index-brem_index]
   excitation_cs = excitation_cross_sections[index]
   cutoff_cs = cutoff_cross_sections[index]
   cutoff_cdf = analog_dist.evaluateCDF( energy, 0.9 )
   analog_cs = cutoff_cs/cutoff_cdf
   ionization_cs = 0.0

   for shell in subshells:
       ionization_cross_section = native_data.getElectroionizationCrossSection(shell)
       ionization_index = native_data.getElectroionizationCrossSectionThresholdEnergyIndex(shell)

       i = index-ionization_index
       if i >= 0:
           ionization_cs += ionization_cross_section[i]

   inelastic_cs = brem_cs + excitation_cs + ionization_cs
   total_cs = inelastic_cs + analog_cs

   print '\tbrem_cs       = ','%.16e' % brem_cs
   print '\texcitation_cs = ','%.16e' % excitation_cs
   print '\tionization_cs = ','%.16e' % ionization_cs
   print '\t-----------------------------------------'
   print '\tinelastic_cs  = ','%.16e' % inelastic_cs
   print '\telastic_cs    = ','%.16e' % analog_cs
   print '\t-----------------------------------------'
   print '\ttotal cs      = ','%.16e' % total_cs




###
### Electroatom Native Factory/Electroatom Factory Test Check
###
print "\n----- Electroatom Factory Classes -----\n"
print "\n----- H -----\n"
tot_elastic_cross_sections = native_data.getTotalElasticCrossSection()
cutoff_cross_sections = native_data.getCutoffElasticCrossSection()
moment_cross_sections = Collision.createLogLogLogExactMomentPreservingElasticReaction(native_data, 0.9, 1e-15)
subshells = native_data.getSubshells()
shell = subshells[0]
ionization_cross_sections = native_data.getElectroionizationCrossSection(shell)
ionization_index = native_data.getElectroionizationCrossSectionThresholdEnergyIndex(shell)
cutoff_dist = Collision.createLinLinLogExactCutoffElasticDistribution(native_data, 1.0, 1e-15)
#print energy_grid
#print moment_cross_sections
energy = 20.0
index = 0
for i in range(0, energy_grid.size ):
   if energy_grid[i] <= energy:
       index = i
print energy_grid[index]
print energy_grid[index+1]
print moment_cross_sections.getCrossSection( energy_grid[index] )
print moment_cross_sections.getCrossSection( energy_grid[index+1] )

energy = 1e-5
brem_cs = brem_cross_sections[0]
excitation_cs = excitation_cross_sections[0]
ionization_cs = ionization_cross_sections[0]
tot_elastic_cs = tot_elastic_cross_sections[0]
tot_cs = brem_cs + excitation_cs + ionization_cs + tot_elastic_cs

cutoff = cutoff_dist.evaluateCDF( energy, 0.9 )
hybrid_cs = cutoff_cross_sections[0]*cutoff + moment_cross_sections.getCrossSection( energy_grid[0] )
tot_cutoff_cs = brem_cs + excitation_cs + ionization_cs + hybrid_cs

max_ionization = ionization_cs/tot_cs
max_elastic = (ionization_cs+tot_elastic_cs)/tot_cs
max_excitation = (ionization_cs+tot_elastic_cs+excitation_cs)/tot_cs
max_brem = (ionization_cs+tot_elastic_cs+excitation_cs+brem_cs)/tot_cs

print "energy = ", energy
print '\tbrem_cs = ','%.16e' % brem_cs
print '\texcitation_cs = ','%.16e' % excitation_cs
print '\tionization_cs = ','%.16e' % ionization_cs
print '\n\ttot_elastic_cs = ','%.16e' % tot_elastic_cs
print '\ttot_cs = ','%.16e' % tot_cs
print '\n\thybrid_cs = ','%.16e' % hybrid_cs
print '\ttot_cutoff_cs = ','%.16e' % tot_cutoff_cs
print '\n\tindex = ', 0
print '\n\tmax ionization random number = ','%.16e' % max_ionization
print '\tmax elastic random number = ','%.16e' % max_elastic
print '\tmax excitation random number = ','%.16e' % max_excitation
print '\tmax brem random number = ','%.16e' % max_brem


energy = 1e-3
index = 0
for i in range(0, energy_grid.size ):
   if energy_grid[i] <= energy:
       index = i

energy_0 = energy_grid[index]
brem_cs_0 = brem_cross_sections[index]
excitation_cs_0 = excitation_cross_sections[index]
ionization_cs_0 = 0.0
if index-ionization_index > 0:
    ionization_cs_0 = ionization_cross_sections[index-ionization_index]
tot_elastic_cs_0 = tot_elastic_cross_sections[index]
cs_0 = brem_cs_0 + excitation_cs_0 + ionization_cs_0 + tot_elastic_cs_0

cutoff_0 = cutoff_dist.evaluateCDF( energy_0, 0.9 )
cutoff_cs_0 = cutoff_cross_sections[index]*cutoff_0
moment_cs_0 = moment_cross_sections.getCrossSection( energy_grid[index] )


energy_1 = energy_grid[index+1]
brem_cs_1 = brem_cross_sections[index+1]
excitation_cs_1 = excitation_cross_sections[index+1]
ionization_cs_1 = 0.0
if index+1-ionization_index > 0:
    ionization_cs_1 = ionization_cross_sections[index+1-ionization_index]

tot_elastic_cs_1 = tot_elastic_cross_sections[index+1]
cs_1 = brem_cs_1 + excitation_cs_1 + ionization_cs_1 + tot_elastic_cs_1

cutoff_1 = cutoff_dist.evaluateCDF( energy_1, 0.9 )
cutoff_cs_1 = cutoff_cross_sections[index+1]*cutoff_1
moment_cs_1 = moment_cross_sections.getCrossSection( energy_grid[index+1] )

brem_cs = brem_cs_0 + (brem_cs_1 - brem_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
excitation_cs = excitation_cs_0 + (excitation_cs_1 - excitation_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
ionization_cs = ionization_cs_0 + (ionization_cs_1 - ionization_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
tot_elastic_cs = tot_elastic_cs_0 + (tot_elastic_cs_1 - tot_elastic_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
tot_cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )

cutoff_cs = cutoff_cs_0 + (cutoff_cs_1 - cutoff_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
moment_cs = moment_cs_0 + (moment_cs_1 - moment_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
hybrid_cs = cutoff_cs+moment_cs
tot_cutoff_cs = brem_cs + excitation_cs + ionization_cs + hybrid_cs

max_excitation = (excitation_cs)/tot_cs
max_brem = (brem_cs+excitation_cs)/tot_cs
max_ionization = (ionization_cs+brem_cs+excitation_cs)/tot_cs
max_elastic = (ionization_cs+brem_cs+excitation_cs+tot_elastic_cs)/tot_cs


print "\nenergy = ", energy
print '\tbrem_cs = ', brem_cs/tot_cs
print '\texcitation_cs = ', excitation_cs/tot_cs
print '\tionization_cs = ', ionization_cs/tot_cs

print '\n\ttot_elastic_cs = ', tot_elastic_cs/tot_cs
print '\ttot_cs = ', tot_cs/tot_cs

print '\n\thybrid_cs = ',hybrid_cs/tot_cs
print '\ttot_cutoff_cs = ', tot_cutoff_cs/tot_cs

print '\n\tindex = ', index

print '\n\tmax excitation random number = ','%.16e' % max_excitation
print '\tmax brem random number = ','%.16e' % max_brem
print '\tmax ionization random number = ','%.16e' % max_ionization
print '\tmax elastic random number = ','%.16e' % max_elastic

#energy = 20.0
#brem_cs = brem_cross_sections[brem_cs.size -1]
#excitation_cs = excitation_cross_sections[excitation_cs.size -1]
#ionization_cs = ionization_cross_sections[ionization_cs.size -1]
#tot_elastic_cs = tot_elastic_cross_sections[tot_elastic_cs.size -1]
#tot_cs = brem_cs + excitation_cs + ionization_cs + tot_elastic_cs

#cutoff = reduced_cutoff_ratio[reduced_cutoff_ratio.size -1]
##print cutoff_dist.evaluateCDF( energy, 0.9 )
##print cutoff

#hybrid_cs = cutoff_cross_sections[cutoff_cs.size -1]*cutoff + moment_cross_sections.getCrossSection( energy_grid[moment_cs.size -1] )
#tot_cutoff_cs = brem_cs + excitation_cs + ionization_cs + hybrid_cs

#max_ionization = ionization_cs/tot_cs
#max_elastic = (ionization_cs+tot_elastic_cs)/tot_cs
#max_excitation = (ionization_cs+tot_elastic_cs+excitation_cs)/tot_cs
#max_brem = (ionization_cs+tot_elastic_cs+excitation_cs+brem_cs)/tot_cs

#print "\nenergy = ", energy
#print '\tbrem_cs = ','%.16e' % brem_cs
#print '\texcitation_cs = ','%.16e' % excitation_cs
#print '\tionization_cs = ','%.16e' % ionization_cs

#print '\n\ttot_elastic_cs = ','%.16e' % tot_elastic_cs
#print '\ttot_cs = ','%.16e' % tot_cs

#print '\n\thybrid_cs = ','%.16e' % hybrid_cs
#print '\ttot_cutoff_cs = ','%.16e' % tot_cutoff_cs

#print '\n\tindex = ', excitation_cs.size -2

#print '\n\tmax ionization random number = ','%.16e' % max_ionization
#print '\tmax elastic random number = ','%.16e' % max_elastic
#print '\tmax excitation random number = ','%.16e' % max_excitation
#print '\tmax brem random number = ','%.16e' % max_brem



####
#### Electroatom Native Factory/Electroatom Factory Test Check
####
#print "\n----- Electroatom Factory Classes -----\n"
#print "\n----- C -----\n"

#data_list = cs_list.get( 'C-Native' )
#file_name = datadir + data_list.get( 'electroatomic_file_path' )
#native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
#energy_grid = native_data.getElectronEnergyGrid()
#brem_cross_sections = native_data.getBremsstrahlungCrossSection()
#excitation_cross_sections = native_data.getAtomicExcitationCrossSection()
#tot_elastic_cross_sections = native_data.getTotalElasticCrossSection()
#cutoff_cross_sections = native_data.getCutoffElasticCrossSection()
#moment_cross_sections = native_data.getMomentPreservingCrossSection()
#subshells = native_data.getSubshells()
#ionization_1_cross_sections = native_data.getElectroionizationCrossSection(subshells[0])
#ionization_2_cross_sections = native_data.getElectroionizationCrossSection(subshells[1])
#ionization_3_cross_sections = native_data.getElectroionizationCrossSection(subshells[2])
#ionization_4_cross_sections = native_data.getElectroionizationCrossSection(subshells[3])

#energy = 1e-5
#brem_cs = brem_cross_sections[0]
#excitation_cs = excitation_cross_sections[0]
#ionization_1_cs = ionization_1_cross_sections[0]
#ionization_2_cs = ionization_2_cross_sections[0]
#ionization_3_cs = ionization_3_cross_sections[0]
#ionization_4_cs = ionization_4_cross_sections[0]
#ionization_cs = ionization_1_cs+ ionization_2_cs + ionization_3_cs + ionization_4_cs
#tot_elastic_cs = tot_elastic_cross_sections[0]
#tot_cs = brem_cs + excitation_cs + ionization_cs + tot_elastic_cs

#cutoff = reduced_cutoff_ratio[0]
#hybrid_cs = cutoff_cross_sections[0]*cutoff + moment_cross_sections[0]
#tot_cutoff_cs = brem_cs + excitation_cs + ionization_1_cs + hybrid_cs

#print "energy = ", energy
#print '\tbrem_cs = ','%.16e' % brem_cs
#print '\texcitation_cs = ','%.16e' % excitation_cs
#print '\tionization_1_cs = ','%.16e' % ionization_1_cs
#print '\tionization_2_cs = ','%.16e' % ionization_2_cs
#print '\tionization_3_cs = ','%.16e' % ionization_3_cs
#print '\tionization_4_cs = ','%.16e' % ionization_4_cs
#print '\n\ttot_elastic_cs = ','%.16e' % tot_elastic_cs
#print '\ttot_cs = ','%.16e' % tot_cs
#print '\n\thybrid_cs = ','%.16e' % hybrid_cs
#print '\ttot_cutoff_cs = ','%.16e' % tot_cutoff_cs
#print '\n\tindex = ', 0

#energy = 1e-3
#index = 0
#for i in range(0, energy_grid.size ):
#    if energy_grid[i] <= energy:
#        index = i

#energy_0 = energy_grid[index]
#brem_cs_0 = brem_cross_sections[index]
#excitation_cs_0 = excitation_cross_sections[index]
#ionization_1_cs_0 = ionization_1_cross_sections[index]
#ionization_2_cs_0 = ionization_2_cross_sections[index]
#ionization_3_cs_0 = ionization_3_cross_sections[index]
#ionization_4_cs_0 = ionization_4_cross_sections[index]
#ionization_cs_0 = ionization_1_cs_0 + ionization_2_cs_0 + ionization_3_cs_0 + ionization_4_cs_0

#tot_elastic_cs_0 = tot_elastic_cross_sections[index]
#cs_0 = brem_cs_0 + excitation_cs_0 + ionization_cs_0 + tot_elastic_cs_0

#cutoff_0 = reduced_cutoff_ratio[index]
#cutoff_cs_0 = cutoff_cross_sections[index]*cutoff_0
#moment_cs_0 = moment_cross_sections[index]


#energy_1 = energy_grid[index+1]
#brem_cs_1 = brem_cross_sections[index+1]
#excitation_cs_1 = excitation_cross_sections[index+1]
#ionization_1_cs_1 = ionization_1_cross_sections[index+1]
#ionization_2_cs_1 = ionization_2_cross_sections[index+1]
#ionization_3_cs_1 = ionization_3_cross_sections[index+1]
#ionization_4_cs_1 = ionization_4_cross_sections[index+1]
#ionization_cs_1 = ionization_1_cs_1 + ionization_2_cs_1 + ionization_3_cs_1 + ionization_4_cs_1
#tot_elastic_cs_1 = tot_elastic_cross_sections[index+1]
#cs_1 = brem_cs_1 + excitation_cs_1 + ionization_cs_1 + tot_elastic_cs_1

#cutoff_1 = reduced_cutoff_ratio[index+1]
#cutoff_cs_1 = cutoff_cross_sections[index+1]*cutoff_1
#moment_cs_1 = moment_cross_sections[index+1]

#brem_cs = brem_cs_0 + (brem_cs_1 - brem_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#excitation_cs = excitation_cs_0 + (excitation_cs_1 - excitation_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#ionization_1_cs = ionization_1_cs_0 + (ionization_1_cs_1 - ionization_1_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#ionization_2_cs = ionization_2_cs_0 + (ionization_2_cs_1 - ionization_2_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#ionization_3_cs = ionization_3_cs_0 + (ionization_3_cs_1 - ionization_3_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#ionization_4_cs = ionization_4_cs_0 + (ionization_4_cs_1 - ionization_4_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#ionization_cs = ionization_1_cs+ ionization_2_cs + ionization_3_cs + ionization_4_cs

#tot_elastic_cs = tot_elastic_cs_0 + (tot_elastic_cs_1 - tot_elastic_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#tot_cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )

#cutoff_cs = cutoff_cs_0 + (cutoff_cs_1 - cutoff_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#moment_cs = moment_cs_0 + (moment_cs_1 - moment_cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
#hybrid_cs = cutoff_cs+moment_cs
#tot_cutoff_cs = brem_cs + excitation_cs + ionization_cs + hybrid_cs

#print "energy = ", energy
#print '\tbrem_cs = ','%.16e' % brem_cs
#print '\texcitation_cs = ','%.16e' % excitation_cs
#print '\tionization_1_cs = ','%.16e' % ionization_1_cs
#print '\tionization_2_cs = ','%.16e' % ionization_2_cs
#print '\tionization_3_cs = ','%.16e' % ionization_3_cs
#print '\tionization_4_cs = ','%.16e' % ionization_4_cs

#print '\n\ttot_elastic_cs = ','%.16e' % tot_elastic_cs
#print '\ttot_cs = ','%.16e' % tot_cs

#print '\n\thybrid_cs = ','%.16e' % hybrid_cs
#print '\ttot_cutoff_cs = ','%.16e' % tot_cutoff_cs

#print '\n\tindex = ', index

#energy = 20.0
#brem_cs = brem_cross_sections[brem_cs.size -1]
#excitation_cs = excitation_cross_sections[excitation_cs.size -1]
#ionization_1_cs = ionization_1_cross_sections[ionization_1_cs.size -1]
#ionization_2_cs = ionization_2_cross_sections[ionization_2_cs.size -1]
#ionization_3_cs = ionization_3_cross_sections[ionization_3_cs.size -1]
#ionization_4_cs = ionization_4_cross_sections[ionization_4_cs.size -1]
#ionization_cs = ionization_1_cs+ ionization_2_cs + ionization_3_cs + ionization_4_cs
#tot_elastic_cs = tot_elastic_cross_sections[tot_elastic_cs.size -1]
#tot_cs = brem_cs + excitation_cs + ionization_cs + tot_elastic_cs

#cutoff = reduced_cutoff_ratio[reduced_cutoff_ratio.size -1]
#hybrid_cs = cutoff_cross_sections[cutoff_cs.size -1]*cutoff + moment_cross_sections[moment_cs.size -1]
#tot_cutoff_cs = brem_cs + excitation_cs + ionization_cs + hybrid_cs

#print "energy = ", energy
#print '\tbrem_cs = ','%.16e' % brem_cs
#print '\texcitation_cs = ','%.16e' % excitation_cs
#print '\tionization_1_cs = ','%.16e' % ionization_1_cs
#print '\tionization_2_cs = ','%.16e' % ionization_2_cs
#print '\tionization_3_cs = ','%.16e' % ionization_3_cs
#print '\tionization_4_cs = ','%.16e' % ionization_4_cs

#print '\n\ttot_elastic_cs = ','%.16e' % tot_elastic_cs
#print '\ttot_cs = ','%.16e' % tot_cs

#print '\n\thybrid_cs = ','%.16e' % hybrid_cs
#print '\ttot_cutoff_cs = ','%.16e' % tot_cutoff_cs

#print '\n\tindex = ', excitation_cs.size -2