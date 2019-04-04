#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.MonteCarlo as MonteCarlo
import PyFrensie.MonteCarlo.Collision as Collision
import PyFrensie.MonteCarlo.Electron as Electron
import numpy
import matplotlib.pyplot as plt

Utility.initFrensiePrng()

# datadir = '/home/software/mcnpdata/native/epr/'
datadir = '/home/lkersting/frensie/src/packages/test_files/native/'

# Set the atomic number (1-99)
atomic_number = 82

# Set the interpolation ("Direct", "Correlated")
interp="Correlated"

# Set the energy (1e-5-1e5)
energy = 1e-2

# energy = 10.5

# Set the xlim and ylim
xlims = [ [0.7,.99], [0.5,1.0] ]
ylims = [ [0.05,0.4], [0.0,1.0] ]

# xlims = [ [0.7,1.0], [0.7,1.0] ]
# ylims = [ [0.0,0.02], [0.0,0.02] ]


### -------------------------------------------------------------------------- ##
###  Hybrid Plots
### -------------------------------------------------------------------------- ##

# filename = datadir + 'epr_native_' + str(atomic_number) + '.xml'
filename = datadir+ 'test_epr_' + str(atomic_number) + '_native.xml'

native_data = Native.ElectronPhotonRelaxationDataContainer( filename )
energy_grid = native_data.getElasticAngularEnergyGrid()

# Distributions
if interp == "Correlated":
  # Distributions
  coupled_dist = Electron.createCoupledElasticDistribution_LogLogCorrelated(native_data, MonteCarlo.TWO_D_UNION, 1e-15)
  cutoff_dist = Electron.createCutoffElasticDistribution_LogLogCorrelated(native_data, 1.0, 1e-15 )
  hybrid_dist = Electron.createHybridElasticDistribution_LogLogCorrelated(native_data, 0.9, 1e-15 )

  # Reactions
  coupled_reaction = Electron.createCoupledElasticReaction_LogLogCorrelated(native_data, MonteCarlo.TWO_D_UNION, 1e-15)
  cutoff_reaction = Electron.createCutoffElasticReaction_LogLogCorrelated(native_data, 1.0, 1e-15 )
  mp_reaction = Electron.createMomentPreservingElasticReaction_LogLogCorrelated(native_data, 0.9, 1e-15 )

elif interp == "Direct":
  # Distributions
  coupled_dist = Electron.createCoupledElasticDistribution_LogLogDirect(native_data, MonteCarlo.TWO_D_UNION, 1e-15)
  cutoff_dist = Electron.createCutoffElasticDistribution_LogLogDirect(native_data, 1.0, 1e-15 )
  hybrid_dist = Electron.createHybridElasticDistribution_LogLogDirect(native_data, 0.9, 1e-15 )

  # Reactions
  coupled_reaction = Electron.createCoupledElasticReaction_LogLogDirect(native_data, MonteCarlo.TWO_D_UNION, 1e-15)
  cutoff_reaction = Electron.createCutoffElasticReaction_LogLogDirect(native_data, 1.0, 1e-15 )
  mp_reaction = Electron.createMomentPreservingElasticReaction_LogLogDirect(native_data, 0.9, 1e-15 )


for i in range(0, len(energy_grid) ):
    if energy_grid[i] < energy:
        index = i
    else:
      break

energy_0 = energy_grid[index]# 8e-3
energy_1 = energy_grid[index+1]# 1.6e-2

cutoff_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )
reduced_cross_section = cutoff_reaction.getCrossSection( energy )*cutoff_cdf
mp_cross_section = mp_reaction.getCrossSection( energy )
tot_cross_section = coupled_reaction.getCrossSection( energy )

print 'energy = ', energy
print 'energy_0 = ', energy_0
print 'energy_0 = ', energy_1
print 'cutoff_cdf = ','%.16e' % cutoff_cdf, "\n"
print 'cutoff cross section  = \t', '%.16e' % cutoff_reaction.getCrossSection( energy )
print "total cross section   = \t", '%.16e' % tot_cross_section
print 'moment preserving cs  = \t', '%.16e' % mp_cross_section
print "reduced_cross_section = \t", '%.16e' % reduced_cross_section, "\n"

cross_section_ratio = reduced_cross_section/mp_cross_section
print 'cross section ratio   = \t','%.16e' % cross_section_ratio
sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
print 'sampling_ratio        = \t','%.16e' % sampling_ratio, "\n"

print "Sample Hybrid at ", energy

upper_random_numbers = numpy.concatenate([[sampling_ratio],numpy.arange(sampling_ratio*1.000001,1.0, 0.0001)])
upper_random_numbers = numpy.append(upper_random_numbers,1.0-1e-15)
below_ratio = sampling_ratio*0.999999

min_cdf_0 = hybrid_dist.evaluateCDF( energy_0, 0.0 )*1.0001
min_cdf_1 = hybrid_dist.evaluateCDF( energy_1, 0.0 )*1.0001
min_cdf = hybrid_dist.evaluateCDF( energy, 0.0 )*1.0001
random_numbers = numpy.append(numpy.arange(min_cdf,sampling_ratio, 0.00001),below_ratio)
random_numbers = numpy.concatenate([random_numbers,upper_random_numbers])


Prng.RandomNumberGenerator.setFakeStream(random_numbers)
sample = [None] * (len( random_numbers ))
for i in range(0,len(random_numbers)):
    sampled_energy, sample[i] = hybrid_dist.sample( energy )

min_cdf_0 = hybrid_dist.evaluateCDF( energy_0, 0.0 )*1.0001
random_numbers_0 = numpy.append(numpy.arange(min_cdf_0,sampling_ratio, 0.00001),below_ratio)
random_numbers_0 = numpy.concatenate([random_numbers_0,upper_random_numbers])

Prng.RandomNumberGenerator.setFakeStream(random_numbers_0)
sample_0 = [None] * (len( random_numbers_0 ))
for i in range(0,len(random_numbers_0)):
    sampled_energy, sample_0[i] = hybrid_dist.sample( energy_0 )

min_cdf_1 = hybrid_dist.evaluateCDF( energy_1, 0.0 )*1.0001
random_numbers_1 = numpy.append(numpy.arange(min_cdf_1,sampling_ratio, 0.00001),below_ratio)
random_numbers_1 = numpy.concatenate([random_numbers_1,upper_random_numbers])

Prng.RandomNumberGenerator.setFakeStream(random_numbers_1)
sample_1 = [None] * (len( random_numbers_1 ))
for i in range(0,len(random_numbers_1)):
    sampled_energy, sample_1[i] = hybrid_dist.sample( energy_1 )

min_cdf_full = coupled_dist.evaluateCDF( energy, 0.0 )*1.0001
random_numbers_full = numpy.append(numpy.arange(min_cdf_full,sampling_ratio, 0.00001),below_ratio)
random_numbers_full = numpy.concatenate([random_numbers_full,upper_random_numbers])

Prng.RandomNumberGenerator.setFakeStream(random_numbers_full)
sample_full = [None] * (len( random_numbers_full ))
for i in range(0,len(random_numbers_full)):
    sampled_energy, sample_full[i] = coupled_dist.sample( energy )


linestyles = [(0, ()), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 5)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

# fig1 = plt.figure(num=1, figsize=(10,5))
# plt.xlabel(r'Angle Cosine ($\mu$)')
# plt.xlabel('CDF')
# plt.title('Hybrid Elastic Scattering Angle vs CDF')
# plt.plot( random_numbers, sample, label=r'$E = 1.0e^{-2}$')
# plt.plot( random_numbers_0, sample_0, label=r'$E_{0} = 8.0e^{-3}$')
# plt.plot( random_numbers_1, sample_1, label=r'$E_{1} = 1.6e^{-2}$')
# # plt.yscale('log')
# plt.ylim(0.85,.99)
# plt.xlim(0.1,0.4)
# plt.legend(loc="best")

# fig1.savefig('../../Desktop/hybrid_cdf_loglog_correlated.pdf', bbox_inches='tight')

# fig2 = plt.figure(num=2, figsize=(10,5))
# plt.xlabel(r'Angle Cosine ($\mu$)')
# plt.xlabel('CDF')
# plt.title('Hybrid Elastic Scattering Angle vs CDF')
# plt.plot( random_numbers, sample, label=r'$E = 1.0e^{-2}$')
# plt.plot( random_numbers_0, sample_0, label=r'$E_{0} = 8.0e^{-3}$')
# plt.plot( random_numbers_1, sample_1, label=r'$E_{1} = 1.6e^{-2}$')
# # plt.yscale('log')
# plt.ylim(0.0,.99)
# plt.xlim(0.0,1.0)
# plt.legend(loc="best")

# fig2.savefig('../../Desktop/hybrid_cdf_loglog_correlated_full.pdf', bbox_inches='tight')

xlim = xlims[0]
ylim = ylims[0]

fig3 = plt.figure(num=3, figsize=(10,5))
plt.xlabel(r'Angle Cosine ($\mu$)')
plt.ylabel('CDF')
plt.title('Hybrid Elastic Scattering Distributions')
plt.plot( sample, random_numbers, color="black", label=r'$E = 1.0e^{-2}$', linewidth=1.8)
plt.plot( sample_0, random_numbers_0, color="red", linestyle="--", label=r'$E_{0} = 8.0e^{-3}$', linewidth=1.8)
plt.plot( sample_1, random_numbers_1, color="green", linestyle="-.", label=r'$E_{1} = 1.6e^{-2}$', linewidth=1.8)
# plt.xscale('log')
plt.xlim(xlim[0], xlim[1])
plt.ylim(ylim[0], ylim[1])
plt.legend(loc=2)

fig3.savefig('./Pb_hybrid_cdf_correlated.pdf', bbox_inches='tight', dpi=600)
fig3.savefig('./Pb_hybrid_cdf_correlated.png', bbox_inches='tight', dpi=600)
fig3.savefig('./Pb_hybrid_cdf_correlated.eps', bbox_inches='tight', dpi=600)

# fig4 = plt.figure(num=4, figsize=(10,5))
# plt.xlabel(r'Angle Cosine ($\mu$)')
# plt.ylabel('CDF')
# plt.title('Hybrid Elastic Scattering Angle vs CDF')
# plt.plot( sample, random_numbers, label=r'$E = 1.0e^{-2}$')
# plt.plot( sample_0, random_numbers_0, label=r'$E_{0} = 8.0e^{-3}$')
# plt.plot( sample_1, random_numbers_1, label=r'$E_{1} = 1.6e^{-2}$')
# # plt.xscale('log')
# plt.xlim(0.5,.99)
# plt.ylim(0.0,0.4)
# plt.legend(loc="best")

# fig4.savefig('../../Desktop/hybrid_cdf_loglog_correlated_switched_full.pdf', bbox_inches='tight')

xlim = xlims[1]
ylim = ylims[1]

hybrid_label = 'Hybrid'
coupled_label = 'Coupled'
fig5 = plt.figure(num=5, figsize=(10,5))
plt.xlabel(r'Angle Cosine ($\mu$)')
plt.ylabel('CDF')
plt.title('Elastic Scattering Distributions at ' + str(energy) + ' MeV')
plt.plot( sample, random_numbers, color="black", label=hybrid_label, linewidth=1.8)
plt.plot( sample_full, random_numbers_full, linestyle="--", label=coupled_label, linewidth=1.8)
plt.legend(loc=2)
# plt.xscale('log')
plt.xlim(xlim[0], xlim[1])
plt.ylim(ylim[0], ylim[1])

fig5.savefig('./Pb_elastic_cdfs_correlated.pdf', bbox_inches='tight', dpi=600)
fig5.savefig('./Pb_elastic_cdfs_correlated.png', bbox_inches='tight', dpi=600)
fig5.savefig('./Pb_elastic_cdfs_correlated.eps', bbox_inches='tight', dpi=600)

plt.show()