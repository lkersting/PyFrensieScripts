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
h_data_list = cs_list.get( 'Pb-Native' )

### -------------------------------------------------------------------------- ##
###  Forward Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
h_native_file_name = datadir + h_data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( h_native_file_name )
energy_grid = native_data.getElectronEnergyGrid()

tot_elastic_cs = native_data.getTotalElasticCrossSection()
cutoff_cs = native_data.getCutoffElasticCrossSection()
screen_rutherford_cs = native_data.getScreenedRutherfordElasticCrossSection()
screen_rutherford_index = native_data.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
moment_cs = native_data.getMomentPreservingCrossSection()
moment_index = native_data.getMomentPreservingCrossSectionThresholdEnergyIndex()
#reduced_cutoff_ratio = native_data.getReducedCutoffCrossSectionRatios()

###
###  Cutoff Distribution/Reaction Unit Test Check
###
#energy = 4e-4
#cutoff_dist = Collision.createAnalogElasticDistribution(h_native_data, False, True, 1e-15)
#analog_dist = Collision.createAnalogElasticDistribution(h_native_data, True, True, 1e-15)
#cutoff_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )
#print 'cutoff_cdf      = ','%.16e' % cutoff_cdf
#cutoff_cdf = analog_dist.evaluateCDF( energy, 0.9 )
#print 'analog_cdf      = ','%.16e' % cutoff_cdf

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy:
#        print 'cutoff_cs_ratio = ','%.16e' %  reduced_cutoff_ratio[i]
#        print 'energy = ', energy_grid[i]
#        print 'cutoff_cs = ','%.16e' %  cutoff_cs[i]
#        print 'moment_preserving_cs = ','%.16e' % moment_cs[i-moment_index]

#print '%.16e' % analog_dist.evaluateScreenedRutherfordPDF( 1.0e-4, 1.0, 2.68213671998009)
#print '%.16e' % analog_dist.evaluateScreenedRutherfordPDF( 5.5e-4, 1.0, 2.68213671998009)
#print '%.16e' % cutoff_dist.evaluateScreenedRutherfordPDF( 5.5e-4, 1.0, 2.68213671998009)
#print '%.16e' % analog_dist.evaluatePDF( 5.5e-4, 1.0)
#print '%.16e' % cutoff_dist.evaluatePDF( 5.5e-4, 1.0)


###
###  Analog Distribution/Reaction Unit Test Check
###
#energy_0 = 8e-3
#energy_1 = 1.6e-2
#energy = 1e-2
#analog_dist = Collision.createAnalogElasticDistribution(native_data, True, True, 1e-15)
#cutoff_cdf = analog_dist.evaluateCDF( energy, 0.9 )
#print 'cutoff_cdf = ','%.16e' % cutoff_cdf
##print energy_grid
#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy:
#        print 'energy = ', energy_grid[i]
#        print energy_grid[i+1]
#        print 'cutoff_cs = ', cutoff_cs[i]
#        print 'moment_preserving_cs = ','%.16e' % moment_cs[i-moment_index]
#        cross_section_ratio = cutoff_cs[i]*cutoff_cdf/moment_cs[i-moment_index]
#        print 'cross section ratio = ','%.16e' % cross_section_ratio
#        sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#        print 'sampling_ratio = ','%.16e' % sampling_ratio

####
####  Hybrid Distribution/Reaction Unit Test Check
####
cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 1.0, "LinLinLog", True, 1e-15)
elastic_energy_grid = native_data.getElasticAngularEnergyGrid()


#angles = native_data.getCutoffElasticAngles(1e-3)
#pdfs = native_data.getCutoffElasticPDF(1e-3)
#print angles[37]
#print pdfs[37]
#print angles[50]
#print pdfs[50]

#energy = 1e-3
#pdf = cutoff_dist.evaluate( energy, 0.0 )
#print 'pdf[0.0] = ','%.16e' % pdf

#pdf = cutoff_dist.evaluate( energy, 0.9 )
#print 'pdf[0.9] = ','%.16e' % pdf,'\n'


#angles = native_data.getCutoffElasticAngles(1e5)
#pdfs = native_data.getCutoffElasticPDF(1e5)
#print angles[7:9]
#print pdfs[7:9]
#print angles[19:21]
#print pdfs[19:21]

#energy = 1e5
#pdf = cutoff_dist.evaluate( energy, 0.0 )
#print '\npdf[0.0] = ','%.16e' % pdf

#pdf = cutoff_dist.evaluate( energy, 0.9 )
#print 'pdf[0.9] = ','%.16e' % pdf

#energy = 1e-5
#cutoff_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )
#print 'cutoff_cdf[','%.6e' % energy,'] = ','%.16e' % cutoff_cdf

#energy_0 = 0.000925526
#cutoff_cdf_0 = cutoff_dist.evaluateCDF( energy_0, 0.9 )
#print 'cutoff_cdf[','%.6e' % energy_0,'] = ','%.16e' % cutoff_cdf_0

#energy = 1e-3
#cutoff_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )
#print 'cutoff_cdf[','%.6e' % energy,'] = ','%.16e' % cutoff_cdf

#energy_1 = 0.00100182
#cutoff_cdf_1 = cutoff_dist.evaluateCDF( energy_1, 0.9 )
#print 'cutoff_cdf[','%.6e' % energy_1,'] = ','%.16e' % cutoff_cdf_1

#energy = 2e-3
#cutoff_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )
#print 'cutoff_cdf[','%.6e' % energy,'] = ','%.16e' % cutoff_cdf

#cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, "LinLinLog", True, 1e-7)
#random_numbers = [ cutoff_cdf_0, cutoff_cdf, cutoff_cdf_1 ]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print cutoff_dist.sample( energy_0 )
#print cutoff_dist.sample( energy )
#print cutoff_dist.sample( energy_1 )

#cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, "LinLinLog", True, 1e-7)

#energy = 1e-4
#sample_ratio = 0.92652534676358189181
#cutoff_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )

#print "Evaluate"
#unorm_evaluate = cutoff_dist.evaluate( energy, 0.9 )
#evaluate = unorm_evaluate*sample_ratio/cutoff_cdf
#print '%.18e' % evaluate

#print "Evaluate PDF"
#unorm_evaluate_pdf = cutoff_dist.evaluatePDF( energy, 0.9 )
#evaluate_pdf = unorm_evaluate_pdf*sample_ratio/cutoff_cdf
#print '%.18e' % evaluate_pdf

#print "Evaluate CDF"
#unorm_evaluate_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )
#evaluate_cdf = unorm_evaluate_cdf*sample_ratio/cutoff_cdf
#print '%.18e' % evaluate_cdf

#cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, "LinLinLog", True, 1e-15)

#energy = 1e-4
#sample_ratio = 0.92974996583144009499
#cutoff_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )
#print cutoff_cdf
#print "Evaluate"
#unorm_evaluate = cutoff_dist.evaluate( energy, 0.9 )
#evaluate = unorm_evaluate*sample_ratio/cutoff_cdf
#print '%.18e' % evaluate

#print "Evaluate PDF"
#unorm_evaluate_pdf = cutoff_dist.evaluatePDF( energy, 0.9 )
#evaluate_pdf = unorm_evaluate_pdf*sample_ratio/cutoff_cdf
#print '%.18e' % evaluate_pdf

#print "Evaluate CDF"
#unorm_evaluate_cdf = cutoff_dist.evaluateCDF( energy, 0.9 )
#evaluate_cdf = unorm_evaluate_cdf*sample_ratio/cutoff_cdf
#print '%.18e' % evaluate_cdf

#discrete_angles_1 = native_data.getMomentPreservingElasticDiscreteAngles(1e-3)
#discrete_weights_1 = native_data.getMomentPreservingElasticWeights(1e-3)
#discrete_energy_grid = native_data.getElasticAngularEnergyGrid()
#print discrete_energy_grid
#print "discrete angles 1e-3"
#print "1: ",'%.18e' % discrete_angles_1[0]
##print "2: ",'%.18e' % discrete_angles_1[1]
#print "discrete weights 1e-3"
#print "1: ",'%.18e' % discrete_weights_1[0]
##print "2: ",'%.18e' % discrete_weights_1[1]

#discrete_angles_0 = native_data.getMomentPreservingElasticDiscreteAngles(1e5)
#discrete_weights_0 = native_data.getMomentPreservingElasticWeights(1e5)
#print "discrete angles 1e-5"
#print "1: ",'%.18e' % discrete_angles_0[0]
##print "2: ",'%.18e' % discrete_angles_0[1]
#print "dicrete weights 1e-5"
#print "1: ",'%.18e' % discrete_weights_0[0]
#print "2: ",'%.18e' % discrete_weights_0[1]

#print "discrete angles 1e-4"
#log_interp = numpy.log(1e-4/1e-5)/numpy.log(1e-3/1e-5)
#lin_interp = (1e-4-1e-5)/(1e-3-1e-5)
#interp = lin_interp
#lower_sample = discrete_angles_0[0]

#angle_1 = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*interp
#angle_2 = discrete_angles_0[0] + (discrete_angles_1[1]-discrete_angles_0[0])*interp
#angle_22 = discrete_angles_1[0] + (discrete_angles_1[1]-discrete_angles_0[0])*interp
#angle_3 = discrete_angles_0[1] + (discrete_angles_1[1]-discrete_angles_0[1])*interp
#print "1: ",'%.18e' % angle_1
#print "2: ",'%.18e' % angle_2
#print "3: ",'%.18e' % angle_3


#print "CrossSectionRatio: ",'%.18e' % hybrid_dist.getCrossSectionRatio(1e-4)
#print "SamplingRatio: ",'%.18e' % hybrid_dist.getSamplingRatio( 1e-4)

#print "Sample Hybrid at ", energy
#random_numbers = numpy.append(numpy.arange(0.0,sampling_ratio, 0.001),sampling_ratio*0.99999)
#random_numbers = numpy.append(random_numbers,sampling_ratio)
#random_numbers = numpy.concatenate([random_numbers,numpy.arange(sampling_ratio*1.001,1.0, 0.001)])
#random_numbers = numpy.append(random_numbers,1.0-1e-15)

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "New Sample Impl"
#new_sample = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    new_sample[i] = hybrid_dist.newSampleImpl( energy )

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "New Sample Impl 2"
#new_sample2 = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    new_sample2[i] = hybrid_dist.newSampleImpl2( energy )

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "Old Sample Impl"
#old_sample = [None] * (len( random_numbers ))
#diff = [None] * (len( random_numbers ))
#diff2 = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    old_sample[i] = hybrid_dist.oldSampleImpl( energy )
#    diff[i] = (old_sample[i] - new_sample[i])#/abs(old_sample[i])
#    diff2[i] = (old_sample[i] - new_sample2[i])#/abs(old_sample[i])

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#old_sample_0 = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    old_sample_0[i] = hybrid_dist.oldSampleImpl( energy_0 )

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#old_sample_1 = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    old_sample_1[i] = hybrid_dist.oldSampleImpl( energy_1 )

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#new_sample_0 = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    new_sample_0[i] = hybrid_dist.newSampleImpl( energy_0 )

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#new_sample_1 = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    new_sample_1[i] = hybrid_dist.newSampleImpl( energy_1 )

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#new_sample2_0 = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    new_sample2_0[i] = hybrid_dist.newSampleImpl2( energy_0 )

#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#new_sample2_1 = [None] * (len( random_numbers ))
#for i in range(0,len(random_numbers)):
#    new_sample2_1[i] = hybrid_dist.newSampleImpl2( energy_1 )

#print diff
#print diff2
print len(random_numbers)
print old_sample[len(old_sample)-2]

fig1 = plt.figure(num=1, figsize=(10,5))
plt.subplot2grid((2,6),(0, 0), colspan=5)
plt.xlabel('Angle Cosine')
plt.ylabel('CDF')
plt.title('Hybrid Elastic Scattering Angle vs CDF')
plt.xlim(0.85,1.0)
#plt.ylim(0.02,0.05)
plt.plot( old_sample, random_numbers, label='Old Sample')
plt.plot( new_sample, random_numbers, label='New Sample')
plt.plot( new_sample2, random_numbers, label='New Sample 2')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.subplot2grid((2,6),(1, 0), colspan=5)
plt.xlabel('Angle Cosine')
plt.ylabel('Error')
plt.xlim(0.85,1.0)
plt.plot( old_sample, diff, label='Error')
plt.plot( old_sample, diff2, label='Error 2')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


fig2 = plt.figure(num=2, figsize=(10,5))
plt.subplot2grid((2,6),(0, 0), colspan=5)
plt.ylabel('Angle Cosine')
plt.xlabel('CDF')
plt.title('Hybrid Elastic Scattering Angle vs CDF')
#plt.xlim(0.85,1.0)
plt.ylim(0.85,1.0)
plt.plot( random_numbers, old_sample, label='Old Sample')
plt.plot( random_numbers, new_sample, label='New Sample')
plt.plot( random_numbers, new_sample2, label='New Sample 2')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.subplot2grid((2,6),(1, 0), colspan=5)
plt.ylabel('Angle Cosine')
plt.xlabel('Error')
#plt.xlim(0.85,1.0)
plt.plot( random_numbers, diff, label='Error')
plt.plot( random_numbers, diff2, label='Error 2')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

fig3 = plt.figure(num=3, figsize=(10,5))
plt.subplot2grid((1,6),(0, 0), colspan=5)
plt.xlabel('Angle Cosine')
plt.ylabel('CDF')
plt.title('Hybrid Elastic Scattering Angle vs CDF for Pb')
plt.xlim(0.8,1.0)
plt.plot( old_sample, random_numbers, label=r'$E = 1e^{-2}$')
plt.plot( old_sample_0, random_numbers, label=r'$E_{0} = 8e^{-3}$')
plt.plot( old_sample_1, random_numbers, label=r'$E_{1} = 1.6e^{-2}$')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

fig3.savefig('../Desktop/hybrid_cdf_linlog2.pdf', bbox_inches='tight')

fig4 = plt.figure(num=4, figsize=(10,5))
plt.subplot2grid((1,6),(0, 0), colspan=5)
plt.ylabel('Angle Cosine')
plt.xlabel('CDF')
plt.title('Hybrid Elastic Scattering Angle vs CDF')
plt.ylim(0.84,.99)
plt.xlim(0.1,0.7)
plt.plot( random_numbers, old_sample, label=r'$E = 1.0e^{-2}$')
plt.plot( random_numbers, old_sample_0, label=r'$E_{0} = 8.0e^{-3}$')
plt.plot( random_numbers, old_sample_1, label=r'$E_{1} = 1.6e^{-2}$')
plt.legend(bbox_to_anchor=(0.97, 0.33), loc=1, borderaxespad=0.)

fig4.savefig('../Desktop/hybrid_cdf_linlog_simple_correlated.pdf', bbox_inches='tight')

fig5 = plt.figure(num=5, figsize=(10,5))
plt.subplot2grid((1,6),(0, 0), colspan=5)
plt.ylabel('Angle Cosine')
plt.xlabel('CDF')
plt.title('Hybrid Elastic Scattering Angle vs CDF')
plt.ylim(0.84,.99)
plt.xlim(0.1,0.7)
plt.plot( random_numbers, new_sample, label=r'$E = 1.0e^{-2}$')
plt.plot( random_numbers, new_sample_0, label=r'$E_{0} = 8.0e^{-3}$')
plt.plot( random_numbers, new_sample_1, label=r'$E_{1} = 1.6e^{-2}$')
plt.legend(bbox_to_anchor=(0.97, 0.33), loc=1, borderaxespad=0.)

fig5.savefig('../Desktop/hybrid_cdf_linlog_modified_correlated.pdf', bbox_inches='tight')

fig6 = plt.figure(num=6, figsize=(10,5))
plt.subplot2grid((1,6),(0, 0), colspan=5)
plt.ylabel('Angle Cosine')
plt.xlabel('CDF')
plt.title('Hybrid Elastic Scattering Angle vs CDF')
plt.ylim(0.84,.99)
plt.xlim(0.1,0.7)
plt.plot( random_numbers, new_sample2, label=r'$E = 1.0e^{-2}$')
plt.plot( random_numbers, new_sample2_0, label=r'$E_{0} = 8.0e^{-3}$')
plt.plot( random_numbers, new_sample2_1, label=r'$E_{1} = 1.6e^{-2}$')
plt.legend(bbox_to_anchor=(0.97, 0.33), loc=1, borderaxespad=0.)

fig5.savefig('../Desktop/hybrid_cdf_linlog_histogram_correlated.pdf', bbox_inches='tight')

#random_numbers = [.26, .26, .26, .26, .26, .26]
#print "Old Sample Impl"
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#sample = hybrid_dist.oldSampleImpl( energy )
#sample_0 = hybrid_dist.oldSampleImpl( energy_0 )
#sample_1 = hybrid_dist.oldSampleImpl( energy_1 )
#interp = sample_0 + (sample_1-sample_0)*numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
#interp2 = sample_0 + (sample_1-sample_0)*(energy-energy_0)/(energy_1-energy_0)
#print 'sample old 0 = ','%.16e' % sample_0
#print 'sample old   = ','%.16e' % sample
#print 'sample interp= ','%.16e' % interp
##print 'sample interp= ','%.16e' % interp2
#print 'sample old 1 = ','%.16e' % sample_1
#print "New Sample Impl"
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#sample = hybrid_dist.newSampleImpl( energy )
#sample_0 = hybrid_dist.newSampleImpl( energy_0 )
#sample_1 = hybrid_dist.newSampleImpl( energy_1 )
#discrete_angles_0 = native_data.getMomentPreservingElasticDiscreteAngles(energy_0)
#discrete_angles_1 = native_data.getMomentPreservingElasticDiscreteAngles(energy_1)

#interp = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*numpy.log10(energy/energy_0)/numpy.log10(energy_1/energy_0)
#interp2 = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*(energy-energy_0)/(energy_1-energy_0)
#print 'sample new 0 = ','%.16e' % sample_0
#print 'sample new   = ','%.16e' % sample
##print 'sample interp= ','%.16e' % interp
##print 'sample interp= ','%.16e' % interp2
#print 'sample new 1 = ','%.16e' % sample_1
#print "New Sample Impl 2"
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#sample = hybrid_dist.newSampleImpl2( energy )
#sample_0 = hybrid_dist.newSampleImpl2( energy_0 )
#sample_1 = hybrid_dist.newSampleImpl2( energy_1 )
#discrete_angles_0 = native_data.getMomentPreservingElasticDiscreteAngles(energy_0)
#discrete_angles_1 = native_data.getMomentPreservingElasticDiscreteAngles(energy_1)

#interp = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*numpy.log10(energy/energy_0)/numpy.log10(energy_1/energy_0)
#interp2 = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*(energy-energy_0)/(energy_1-energy_0)
#print 'sample new 0 = ','%.16e' % sample_0
#print 'sample new   = ','%.16e' % sample
#print 'sample interp= ','%.16e' % interp
##print 'sample interp= ','%.16e' % interp2
#print 'sample new 1 = ','%.16e' % sample_1

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy_0:
#        print "\n",energy_0
#        cross_section_ratio = cutoff_cs[i]*cutoff_cdf/moment_cs[i-moment_index]
#        print 'cross section ratio = ','%.16e' % cross_section_ratio
#        sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#        print 'sampling_ratio = ','%.16e' % sampling_ratio
#        print 'scaled random number 1 = ','%.16e' % float((1.0+cross_section_ratio)*0.26)
#        print 'scaled random number 2 = ','%.16e' % float(0.26*(1.0 + cross_section_ratio) - cross_section_ratio)

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy:
#        print "\n",energy
#        cross_section_ratio = cutoff_cs[i]*cutoff_cdf/moment_cs[i-moment_index]
#        print 'cross section ratio = ','%.16e' % cross_section_ratio
#        sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#        print 'sampling_ratio = ','%.16e' % sampling_ratio
#        print 'scaled random number 1 = ','%.16e' % float((1.0+cross_section_ratio)*0.26)
#        print 'scaled random number 2 = ','%.16e' % float(0.26*(1.0 + cross_section_ratio) - cross_section_ratio)

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy_1:
#        print "\n",energy_1
#        cross_section_ratio = cutoff_cs[i]*cutoff_cdf/moment_cs[i-moment_index]
#        print 'cross section ratio = ','%.16e' % cross_section_ratio
#        sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#        print 'sampling_ratio = ','%.16e' % sampling_ratio
#        print 'scaled random number 1 = ','%.16e' % float((1.0+cross_section_ratio)*0.26)
#        print 'scaled random number 2 = ','%.16e' % float(0.26*(1.0 + cross_section_ratio) - cross_section_ratio)

#print "discrete_angles = ", native_data.getMomentPreservingElasticDiscreteAngles(energy_0)
#print "discrete_weights = ",native_data.getMomentPreservingElasticWeights(energy_0)

#print "discrete_angles = ", native_data.getMomentPreservingElasticDiscreteAngles(energy_1)
#print "discrete_weights = ",native_data.getMomentPreservingElasticWeights(energy_1)

#plt.show()

#discrete_angles = h_native_data.getMomentPreservingElasticDiscreteAngles(energy_0)
#discrete_weights = h_native_data.getMomentPreservingElasticWeights(energy_0)
#discrete_energy_grid = h_native_data.getElasticAngularEnergyGrid()
#print discrete_angles
#print discrete_weights

## -------------------------------------------------------------------------- ##
##  Adjoint Elastic Unit Test Check
## -------------------------------------------------------------------------- ##
#h_adjoint_file_name = datadir + h_data_list.get( 'adjoint_electroatomic_file_path' )
#h_adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( h_adjoint_file_name )
#adjoint_energy_grid = h_adjoint_data.getAdjointElectronEnergyGrid()


#tot_adjoint_elastic_cs = h_adjoint_data.getAdjointTotalElasticCrossSection()
#adjoint_cutoff_cs = h_adjoint_data.getAdjointCutoffElasticCrossSection()
#reduced_cutoff_ratio = h_adjoint_data.getReducedCutoffCrossSectionRatios()
#adjoint_screen_rutherford_cs = h_adjoint_data.getAdjointScreenedRutherfordElasticCrossSection()
#adjoint_screen_rutherford_index = h_adjoint_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
#adjoint_moment_cs = h_adjoint_data.getAdjointMomentPreservingCrossSection()
#adjoint_moment_index = h_adjoint_data.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex()

#adjoint_excitation_cs = h_adjoint_data.getAdjointAtomicExcitationCrossSection()
#adjoint_ionization_cs = h_adjoint_data.getAdjointElectroionizationCrossSection(1)
#adjoint_brem_cs = h_adjoint_data.getAdjointBremsstrahlungElectronCrossSection()


#discrete_angles = h_adjoint_data.getAdjointMomentPreservingElasticDiscreteAngles(1e-3)
#discrete_weights = h_adjoint_data.getAdjointMomentPreservingElasticWeights(1e-3)
#discrete_energy_grid = native_data.getElasticAngularEnergyGrid()
#print discrete_energy_grid
#print "discrete angles 1e-3"
#print "1: ",'%.18e' % discrete_angles[0]
#print "2: ",'%.18e' % discrete_angles[1]
#print "discrete weights 1e-3"
#print "1: ",'%.18e' % discrete_weights[0]
#print "2: ",'%.18e' % discrete_weights[1]

#discrete_angles = h_adjoint_data.getAdjointMomentPreservingElasticDiscreteAngles(20.0)
#discrete_weights = h_adjoint_data.getAdjointMomentPreservingElasticWeights(20.0)
#discrete_energy_grid = native_data.getElasticAngularEnergyGrid()
#print discrete_energy_grid
#print "discrete angles 20.0"
#print "1: ",'%.18e' % discrete_angles[0]
#print "2: ",'%.18e' % discrete_angles[1]
#print "discrete weights 20.0"
#print "1: ",'%.18e' % discrete_weights[0]
#print "2: ",'%.18e' % discrete_weights[1]


###
###  Cutoff Distribution/Reaction Unit Test Check
###
#adjoint_cutoff_dist = Collision.createCutoffElasticDistribution(h_adjoint_data, 1.0, True, True, 1e-7)
#energy = 1e-5
#cutoff_cdf = adjoint_cutoff_dist.evaluateCDF( energy, 0.9 )
#print '\n\ncutoff_cdf = ','%.16e' % cutoff_cdf

#energy = 1e-3
#cutoff_cdf = adjoint_cutoff_dist.evaluateCDF( energy, 0.9 )
#print 'cutoff_cdf = ','%.16e' % cutoff_cdf

#energy = 1e5
#cutoff_cdf = adjoint_cutoff_dist.evaluateCDF( energy, 0.9 )
#print 'cutoff_cdf = ','%.16e' % cutoff_cdf
#ratio = reduced_cutoff_ratio[reduced_cutoff_ratio.size -1]
#cutoff_cs = adjoint_cutoff_cs[reduced_cutoff_ratio.size -1]
#moment_cs = adjoint_moment_cs[reduced_cutoff_ratio.size -1]
#hybrid_cs = cutoff_cdf*cutoff_cs + moment_cs

#excitation_cs = adjoint_excitation_cs[adjoint_excitation_cs.size -1]
#ionization_cs = adjoint_ionization_cs[adjoint_ionization_cs.size -1]
#brem_cs = adjoint_brem_cs[adjoint_brem_cs.size -1]
#total_cs = hybrid_cs + excitation_cs + ionization_cs + brem_cs

#print 'cutoff_cdf = ','%.16e' % cutoff_cdf
#print 'ratio      = ','%.16e' % ratio
#print 'cutoff_cs = ','%.16e' % cutoff_cs
#print 'moment_cs = ','%.16e' % moment_cs
#print 'hybrid_cs = ','%.16e' % hybrid_cs

#print 'excitation_cs = ','%.16e' % excitation_cs
#print 'ionization_cs = ','%.16e' % ionization_cs
#print 'brem_cs = ','%.16e' % brem_cs

#print 'total_cs = ','%.16e' % total_cs

#print 'excitation_cs(1e-5) = ','%.16e' % adjoint_excitation_cs[0]

#energy = 1e-3
#energy_1 = adjoint_energy_grid[42]
#energy_2 = adjoint_energy_grid[43]
#linlin_slope = ( energy - energy_1 )/( energy_2 - energy_1 )
#excitation_cs_1 = adjoint_excitation_cs[42]
#excitation_cs_2 = adjoint_excitation_cs[43]

#excitation_cs = excitation_cs_1 + ( excitation_cs_2 - excitation_cs_1 )*linlin_slope
#print 'excitation_cs(1e-3) = ','%.16e' % excitation_cs

###
###  Analog Distribution/Reaction Unit Test Check
###
#adjoint_analog_dist = Collision.createAnalogElasticDistribution(h_adjoint_data)
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

###
###  Hybrid Distribution/Reaction Unit Test Check
###
#adjoint_analog_dist = Collision.createHybridElasticDistribution(h_adjoint_data, 0.9)
#print "Sample Adjoint Hybrid"
#random_numbers = [2.8368270631620132e-02, 1.5494337020438850e-01, 1.5495e-01, 2.9223E-01]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print adjoint_analog_dist.sample( 1.0e-3 )
#print adjoint_analog_dist.sample( 1.0e-3 )
#print adjoint_analog_dist.sample( 1.0e-3 )
#print adjoint_analog_dist.sample( 1.0e-3 )
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

