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
atomic_number = 13

# Set the interpolation ("Direct", "Correlated")
interp="Correlated"

# Set the energy (1e-5-1e5)
energy = 1e-2

energy = 15.7

# Set the xlim and ylim
xlims = [ [0.7,.99], [0.5,1.0] ]
ylims = [ [0.05,0.4], [0.0,1.0] ]

xlims = [ [0.7,1.0], [0.7,1.0] ]
ylims = [ [0.0,0.02], [0.0,0.02] ]


### -------------------------------------------------------------------------- ##
###  Coupled Plots
### -------------------------------------------------------------------------- ##

# filename = datadir + 'epr_native_' + str(atomic_number) + '.xml'
filename = datadir+ 'test_epr_' + str(atomic_number) + '_native.xml'

native_data = Native.ElectronPhotonRelaxationDataContainer( filename )
energy_grid = native_data.getElasticAngularEnergyGrid()

# Distributions
if interp == "Correlated":
  # Distributions
  coupled_dist1 = Electron.createCoupledElasticDistribution_LogLogCorrelated(native_data, MonteCarlo.ONE_D_UNION, 1e-15)
  coupled_dist2 = Electron.createCoupledElasticDistribution_LogLogCorrelated(native_data, MonteCarlo.TWO_D_UNION, 1e-15)
  coupled_dist3 = Electron.createCoupledElasticDistribution_LogLogCorrelated(native_data, MonteCarlo.MODIFIED_TWO_D_UNION, 1e-15)

  # Reactions
  coupled_reaction = Electron.createCoupledElasticReaction_LogLogCorrelated(native_data, MonteCarlo.TWO_D_UNION, 1e-15)

elif interp == "Direct":
  # Distributions
  coupled_dist = Electron.createCoupledElasticDistribution_LogLogDirect(native_data, MonteCarlo.ONE_D_UNION, 1e-15)
  coupled_dist2 = Electron.createCoupledElasticDistribution_LogLogDirect(native_data, MonteCarlo.TWO_D_UNION, 1e-15)
  coupled_dist3 = Electron.createCoupledElasticDistribution_LogLogDirect(native_data, MonteCarlo.MODIFIED_TWO_D_UNION, 1e-15)

  # Reactions
  coupled_reaction = Electron.createCoupledElasticReaction_LogLogDirect(native_data, MonteCarlo.TWO_D_UNION, 1e-15)

print "Sample Coupled at ", energy

min_cdf = coupled_dist1.evaluateCDF( energy, 0.0 )*1.001
step = 1e-5
random_numbers = numpy.append(numpy.arange(0.0,1.0-1e-15, step),1.0-1e-15)

Prng.RandomNumberGenerator.setFakeStream(random_numbers)
sample1 = [None] * (len( random_numbers ))
for i in range(0,len(random_numbers)):
    sampled_energy, sample1[i] = coupled_dist1.sample( energy )

Prng.RandomNumberGenerator.setFakeStream(random_numbers)
sample2 = [None] * (len( random_numbers ))
for i in range(0,len(random_numbers)):
    sampled_energy, sample2[i] = coupled_dist2.sample( energy )

Prng.RandomNumberGenerator.setFakeStream(random_numbers)
sample3 = [None] * (len( random_numbers ))
for i in range(0,len(random_numbers)):
    sampled_energy, sample3[i] = coupled_dist3.sample( energy )


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

max_diff12 = 0.0
max_diff13 = 0.0
max_diff23 = 0.0
index12 = 0
index13 = 0
index23 = 0
for i in range( len(sample1)):
  diff12 = abs(sample1[i] - sample2[i])
  diff13 = abs(sample1[i] - sample3[i])
  diff23 = abs(sample2[i] - sample3[i])
  if diff12 > 1e-8 or diff13 > 1e-8:
    print random_numbers[i], sample1[i], diff12, diff13
  elif diff23 > 1e-8:
    print random_numbers[i], sample2[i], diff23

  if max_diff12 < diff12:
    max_diff12 = diff12
    index12 = i
  if max_diff13 < diff13:
    max_diff13 = diff13
    index13 = i
  if max_diff23 < diff23:
    max_diff23 = diff23
    index23 = i

print random_numbers[index12], sample1[index12], max_diff12
print random_numbers[index13], sample1[index13], max_diff13
print random_numbers[index23], sample2[index23], max_diff23




xlim = xlims[0]
ylim = ylims[0]

fig3 = plt.figure(num=3, figsize=(10,5))
plt.xlabel(r'Angle Cosine ($\mu$)')
plt.ylabel('CDF')
plt.title('Coupled Elastic Scattering Distributions')
plt.plot( sample1, random_numbers, color="black", label='1D Union', linewidth=1.8)
plt.plot( sample2, random_numbers, color="red", linestyle="--", label='2D Union', linewidth=1.8)
plt.plot( sample3, random_numbers, color="green", linestyle="-.", label='M2D Union', linewidth=1.8)
plt.xscale('log')
plt.yscale('log')
# plt.xlim(xlim[0], xlim[1])
# plt.ylim(ylim[0], ylim[1])
plt.ylim(9e-6, 1.0)
plt.legend(loc=2)

output = 'Pb_elastic_cdf_correlated'
fig3.savefig( output + '.pdf', bbox_inches='tight', dpi=600)
fig3.savefig( output + '.png', bbox_inches='tight', dpi=600)
fig3.savefig( output + '.eps', bbox_inches='tight', dpi=600)

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

# xlim = xlims[1]
# ylim = ylims[1]

# hybrid_label = 'Hybrid'
# coupled_label = 'Coupled'
# fig5 = plt.figure(num=5, figsize=(10,5))
# plt.xlabel(r'Angle Cosine ($\mu$)')
# plt.ylabel('CDF')
# plt.title('Elastic Scattering Distributions at ' + str(energy) + ' MeV')
# plt.plot( sample, random_numbers, color="black", label=hybrid_label, linewidth=1.8)
# plt.plot( sample_full, random_numbers_full, linestyle="--", label=coupled_label, linewidth=1.8)
# plt.legend(loc=2)
# # plt.xscale('log')
# plt.xlim(xlim[0], xlim[1])
# plt.ylim(ylim[0], ylim[1])

# fig5.savefig('./Pb_elastic_cdfs_correlated.pdf', bbox_inches='tight', dpi=600)
# fig5.savefig('./Pb_elastic_cdfs_correlated.png', bbox_inches='tight', dpi=600)
# fig5.savefig('./Pb_elastic_cdfs_correlated.eps', bbox_inches='tight', dpi=600)

plt.show()