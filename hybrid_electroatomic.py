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

## -------------------------------------------------------------------------- ##
##  Elastic Data
## -------------------------------------------------------------------------- ##
elements = ['Pb-Native','Al-Native','H-Native']
interps = ["LogLogLog", "LinLinLin", "LinLinLog"]
energies = [1e-5, 1e-3, 1e5 ]

elements = ['Pb-Native']
interps = ["LogLogLog"]
energies = [1e-5, 1.995260e-4, 4e-4, 1e-3, 2e-1, 1e5 ]

for z in elements:
    print "\n----------------------------"
    print "-----", z, "Tests -----"
    print "----------------------------"
    data_list = cs_list.get( z )
    file_name = datadir + data_list.get( 'electroatomic_file_path' )
    native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
    energy_grid = native_data.getElectronEnergyGrid()


    tot_elastic_cs = native_data.getTotalElasticCrossSection()
    cutoff_cross_sections = native_data.getCutoffElasticCrossSection()
    screen_rutherford_cs = native_data.getScreenedRutherfordElasticCrossSection()
    screen_rutherford_index = native_data.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
    moment_cross_sections = Collision.createLogLogLogExactMomentPreservingElasticReaction(native_data, 0.9, 1e-7)


#    if z == 'Pb-Native':
#        energies = [1e-5,4e-4,1e5]
#    else:
#        energies = [1e-5, 1e-3, 1e5 ]

    for interp in interps:
        print "\n--- ",interp,"Tests ---"

        cutoff_dist = Collision.createLogLogLogExactCutoffElasticDistribution(native_data, 0.9, 1e-7)
        if interp == "LinLinLog":
          cutoff_dist = Collision.createLinLinLogExactCutoffElasticDistribution(native_data, 0.9, 1e-15)
        elif interp == "LinLinLin":
          cutoff_dist = Collision.createLinLinLinExactCutoffElasticDistribution(native_data, 0.9, 1e-15)

        ###
        ###  Moment Preserving Reaction Unit Test Check
        ###
        for energy in energies:

            index = 0
            for i in range(0, energy_grid.size ):
                if energy_grid[i] <= energy:
                    index = i

            energy_0 = energy_grid[index]
            moment_cs_0 = moment_cross_sections.getCrossSection( energy_0 )
            cutoff_cdf_0 = cutoff_dist.evaluateCutoffCrossSectionRatio( energy_0 )
            cutoff_cs_0 = cutoff_cross_sections[index]*cutoff_cdf_0
            cs_0 = moment_cs_0 + cutoff_cs_0
            cutoff_cs = cutoff_cs_0
            moment_cs = moment_cs_0
            cs = cs_0


            if energy_0 != energy:
                energy_1 = energy_grid[index+1]
                moment_cs_1 = moment_cross_sections.getCrossSection( energy_1 )
                cutoff_cdf_1 = cutoff_dist.evaluateCutoffCrossSectionRatio( energy_1 )
                cutoff_cs_1 = cutoff_cross_sections[index+1]*cutoff_cdf_1
                cs_1 = moment_cs_1 + cutoff_cs_1

                log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
                lin_interp = (energy-energy_0)/(energy_1-energy_0)

                if interp =="LogLogLog":
                    cutoff_cs = (cutoff_cs_0)*pow((cutoff_cs_1)/(cutoff_cs_0),log_interp)
                    moment_cs = (moment_cs_0)*pow((moment_cs_1)/(moment_cs_0),log_interp)
                    cs = (cs_0)*pow((cs_1)/(cs_0),log_interp)
                elif interp == "LinLinLin":
                    cutoff_cs = cutoff_cs_0 + (cutoff_cs_1-cutoff_cs_0)*lin_interp
                    moment_cs = moment_cs_0 + (moment_cs_1-moment_cs_0)*lin_interp
                    cs = cs_0 + (cs_1-cs_0)*lin_interp
                else:
                    cutoff_cs = cutoff_cs_0 + (cutoff_cs_1-cutoff_cs_0)*log_interp
                    moment_cs = moment_cs_0 + (moment_cs_1-moment_cs_0)*log_interp
                    cs = cs_0 + (cs_1-cs_0)*log_interp
            #moment_cs = moment_cross_sections.getCrossSection( energy )
            print '   Energy = ',energy
            print '\tcutoff cs = ','%.16e' % cutoff_cs
            print '\tmoment cs = ','%.16e' % moment_cs
            print '\thybrid cs = ','%.16e' % cs



#data_list = cs_list.get( 'Al-Native' )
#file_name = datadir + data_list.get( 'electroatomic_file_path' )
#native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
#energy_grid = native_data.getElectronEnergyGrid()


#tot_elastic_cs = native_data.getTotalElasticCrossSection()
#forward_cutoff_cs = native_data.getCutoffElasticCrossSection()
#screen_rutherford_cs = native_data.getScreenedRutherfordElasticCrossSection()
#screen_rutherford_index = native_data.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
#forward_moment_cs = native_data.getMomentPreservingCrossSection()
#forward_moment_index = native_data.getMomentPreservingCrossSectionThresholdEnergyIndex()

####
####  Hybrid Reaction Unit Test Check
####
#cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, "LinLinLin", True, 1e-15)

#energies = [1e-5, 1e-3, 1e5 ]

#for energy in energies:
#    print "\nenergy = ", energy

#    index = 0
#    for i in range(0, energy_grid.size ):
#        if energy_grid[i] <= energy:
#            index = i

#    moment_cs = forward_moment_cross_sections[index]
#    cutoff_cs = forward_cutoff_cross_sections[index]
#    ratio = cutoff_dist.evaluateCutoffCrossSectionRatio( energy )
#    hybrid_cs = cutoff_cs*ratio + moment_cs

#    print "\tdiscrete cross section = " '%.18e' % moment_cs
#    print "\tcutoff cross section   = " '%.18e' % cutoff_cs
#    print "\treduced cutoff ratio   = " '%.18e' % ratio
#    print "\thybrid cross section   = " '%.18e' % hybrid_cs

#energy = 1e-3
#print "\nenergy = ", energy
#index = 0
#for i in range(0, energy_grid.size ):
#    if energy_grid[i] < energy:
#        index = i

#energy_0 = energy_grid[index]
#print "   energy_0 = ", energy_0
#moment_cs_0 = forward_moment_cross_sections[index-forward_moment_index]
#cutoff_cs_0 = forward_cutoff_cross_sections[index]
#ratio_0 = cutoff_dist.evaluateCutoffCrossSectionRatio( energy_0 )
#hybrid_cs_0 = cutoff_cs_0*ratio_0 + moment_cs_0
#print "   hybrid_cs_0 = ", hybrid_cs_0, "\n"

#energy_1 = energy_grid[index+1]
#print "   energy_1 = ", energy_1
#moment_cs_1 = forward_moment_cross_sections[index+1-forward_moment_index]
#cutoff_cs_1 = forward_cutoff_cross_sections[index+1]
#ratio_1 = cutoff_dist.evaluateCutoffCrossSectionRatio( energy_1 )
#hybrid_cs_1 = cutoff_cs_1*ratio_1 + moment_cs_1
#print "   hybrid_cs_1 = ", hybrid_cs_1


#cs = hybrid_cs_0 + ( hybrid_cs_1 - hybrid_cs_0 )*( energy - energy_0 )/( energy_1 - energy_0 )

#print "hybrid cross section = " '%.18e' % cs


#energy = 1e5
#print "\nenergy = ", energy
#moment_cs = forward_moment_cross_sections[forward_moment_cs.size-1]
#cutoff_cs = forward_cutoff_cross_sections[forward_cutoff_cs.size-1]
#ratio = cutoff_dist.evaluateCutoffCrossSectionRatio( energy )
#hrbrid_cs = cutoff_cs*ratio + moment_cs

#print "discrete cross section = " '%.18e' % moment_cs
#print "cutoff cross section = " '%.18e' % cutoff_cs
#print "reduced cutoff ratio = " '%.18e' % ratio
#print "hybrid cross section = " '%.18e' % hrbrid_cs


#analog_dist = Collision.createHybridElasticDistribution(native_data, 0.9, "LogLogLog", True, 1e-15)



#print '%.16e' % analog_dist.evaluateCDF( 1.0e-3, 0.54 )
#print '%.16e' % analog_dist.evaluateCDF( 1.0e-3, 0.9 )
#print '%.16e' % analog_dist.evaluateCDF( 1.0e-3, 0.9239 )
#print '%.16e' % analog_dist.evaluateCDF( 1.0e-3, 0.9788926224755288 )

#discrete_angles = native_data.getMomentPreservingElasticDiscreteAngles(1e-3)
#discrete_weights = native_data.getMomentPreservingElasticWeights(1e-3)
#discrete_energy_grid = native_data.getElasticAngularEnergyGrid()
#print discrete_angles
#print discrete_weights

#energy = 1e-3
#energy_1 = energy_grid[42]
#ratio_1 = cutoff_dist.evaluateCutoffCrossSectionRatio( energy_1 )
#cutoff_1 = forward_cutoff_cross_sections[42]
#moment_1 = forward_moment_cross_sections[42]
#cross_section_1 = cutoff_1*ratio_1 + moment_1

#energy_2 = energy_grid[43]
#ratio_2 = cutoff_dist.evaluateCutoffCrossSectionRatio( energy_2 )
#cutoff_2 = forward_cutoff_cross_sections[43]
#moment_2 = forward_moment_cross_sections[43]
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


#### -------------------------------------------------------------------------- ##
####  Forward Elastic Unit Test Check
#### -------------------------------------------------------------------------- ##
#data_list = cs_list.get( 'Al-Native' )
#native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
#native_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )
#energy_grid = native_data.getElectronEnergyGrid()

#tot_elastic_cs = native_data.getTotalElasticCrossSection()
#cutoff_cs = native_data.getCutoffElasticCrossSection()
#screen_rutherford_cs = native_data.getScreenedRutherfordElasticCrossSection()
#screen_rutherford_index = native_data.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
#moment_cs = native_data.getMomentPreservingCrossSection()
#moment_index = native_data.getMomentPreservingCrossSectionThresholdEnergyIndex()
#discrete_energy_grid = native_data.getElasticAngularEnergyGrid()


####
####  Evaluate Distribution
####
#print '--- Evaluate At Cutoff ---'
#hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, "LinLinLog", True, 1e-7)
#cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, "LinLinLog", True, 1e-7)

#energies = [1e-4, 1e-3, 1e5 ]

#for energy in energies:
#    ratio = hybrid_dist.getSamplingRatio( energy )
#    unorm_eval = cutoff_dist.evaluate(energy, 0.9)
#    evaluation = unorm_eval*ratio

#    print '\nenergy =', energy
#    print '\tsampling ratio =','%.16e' % ratio
#    print '\tunorm eval =','%.16e' % unorm_eval
#    print '\teval   =','%.16e' % evaluation
#    print '\tresult =','%.16e' % hybrid_dist.evaluate( energy, 0.9 )


####
####  Evaluate PDF
####
#print '\n--- Evaluate PDF At Cutoff ---'

#energies = [1e-4, 1e-3, 1e5 ]

#for energy in energies:
#    ratio = hybrid_dist.getSamplingRatio( energy )
#    unorm_pdf = cutoff_dist.evaluatePDF(energy, 0.9)
#    pdf = unorm_pdf*ratio

#    print '\nenergy =', energy
#    print '\tsampling ratio =','%.16e' % ratio
#    print '\tunorm pdf =','%.16e' % unorm_pdf
#    print '\tpdf   =','%.16e' % pdf
#    print '\tresult =','%.16e' % hybrid_dist.evaluatePDF( energy, 0.9 )

####
####  Evaluate CDF
####
#print '\n--- Evaluate CDF At Cutoff ---'

#energies = [1e-4, 1e-3, 1e5 ]

#for energy in energies:
#    ratio = hybrid_dist.getSamplingRatio( energy )
#    unorm_cdf = cutoff_dist.evaluateCDF(energy, 0.9)
#    cdf = unorm_cdf*ratio

#    print '\nenergy =', energy
#    print '\tsampling ratio =','%.16e' % ratio
#    print '\tunorm cdf =','%.16e' % unorm_cdf
#    print '\tcdf   =','%.16e' % cdf
#    print '\tresult =','%.16e' % hybrid_dist.evaluateCDF( energy, 0.9 )


####
####  Cutoff Distribution
####
#print '\n--- Sampling ---'

#energy = 1e-3
#partial_cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, "LinLinLog", True, 1e-14)
#cutoff_cdf = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy )

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy:
#        print '\nenergy = ', energy_grid[i]
#        print '\tcutoff_cdf = ','%.16e' % cutoff_cdf
#        print '\tcutoff_cs = ','%.16e' % cutoff_cross_sections[i]
#        print '\tmoment_preserving_cs = ','%.16e' % moment_cross_sections[i-moment_index]
#        cross_section_ratio = cutoff_cross_sections[i]*cutoff_cdf/moment_cross_sections[i-moment_index]
#        print '\tcross section ratio = ','%.16e' % cross_section_ratio
#        sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#        print '\tsampling_ratio = ','%.16e' % sampling_ratio

#print "\n\t--- continuous sampling ---"
#cdf = partial_cutoff_dist.evaluateCDF( energy, 0.1 )
#print '\tcdf at 0.1 for cutoff = ','%.16e' % cdf
#cdf_scaled = cdf*sampling_ratio
#print '\tcdf at 0.1 for hybrid = ','%.16e' % cdf_scaled
#print '\tcdf at 0.9 for cutoff = ','%.16e' % cutoff_cdf
#print '\tcdf at 0.9 for hybrid = ','%.16e' % sampling_ratio

#discrete_angles = native_data.getMomentPreservingElasticDiscreteAngles(energy)
#discrete_weights = native_data.getMomentPreservingElasticWeights(energy)

#print "\n\t--- discrete sampling ---"
#print '\tdiscrete_weights[0] = ','%.16e' % discrete_weights[0]
#print '\tdiscrete_angles[0] = ','%.16e' % discrete_angles[0]
#weight_0_ratio = (discrete_weights[0] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\tweight_0 cdf range = ','%.16e' % sampling_ratio,' < random number <= ','%.16e' % weight_0_ratio
#print '\tdiscrete_weights[1] = ','%.16e' % discrete_weights[1]
#print '\tdiscrete_angles[1] = ','%.16e' % discrete_angles[1]
#weight_1_ratio = (discrete_weights[0] + discrete_weights[1] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\tweight_1 cdf range = ','%.16e' % weight_0_ratio,' < random number <= ','%.16e' % weight_1_ratio


#energy = 1e-4
#hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, "LinLinLog", True, 1e-7)
#mp_dist = Collision.createMomentPreservingElasticDistribution(native_data, 0.9, "LinLinLog", True, 1e-15)
#cutoff_cdf = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy)

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy:
#        print '\nenergy = ', energy_grid[i]
#        print '\tcutoff_cdf = ','%.16e' % cutoff_cdf
#        print '\tcutoff_cs = ','%.16e' % cutoff_cross_sections[i]
#        print '\tmoment_preserving_cs = ','%.16e' % moment_cross_sections[i-moment_index]
#        cross_section_ratio = cutoff_cross_sections[i]*cutoff_cdf/moment_cross_sections[i-moment_index]
#        print '\tcross section ratio = ','%.16e' % cross_section_ratio
#        cross_section_ratio = hybrid_dist.getCrossSectionRatio( energy )
#        print '\tcross section ratio = ','%.16e' % cross_section_ratio
#        sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#        print '\tsampling_ratio = ','%.16e' % sampling_ratio


#print "\n\t--- continuous sampling ---"
#cdf = partial_cutoff_dist.evaluateCDF( energy, 0.1 )
#print '\tcdf at 0.1 for partial cutoff = ','%.16e' % cdf
#cdf_scaled = cdf*sampling_ratio
#print '\tcdf at 0.1 for hybrid = ','%.16e' % cdf_scaled
#print '\tcdf at 0.9 for cutoff = ','%.16e' % cutoff_cdf
#print '\tcdf at 0.9 for hybrid = ','%.16e' % sampling_ratio

#random_numbers = [cdf,cdf_scaled]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "\tsample[",energy,",",random_numbers[0],"] = ", partial_cutoff_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[1],"] = ",hybrid_dist.sample( energy )

#print "\n\t--- discrete sampling ---"

#index = 0
#for i in range(0, discrete_energy_grid.size ):
#    if discrete_energy_grid[i] <= energy:
#        index = i
#energy_0 = discrete_energy_grid[index]
#energy_1 = discrete_energy_grid[index+1]

#discrete_angles_0 = native_data.getMomentPreservingElasticDiscreteAngles(energy_0)
#discrete_weights_0 = native_data.getMomentPreservingElasticWeights(energy_0)
#print "\n\t--- discrete angles ",energy_0," ---"
#print "\t1: ",'%.18e' % discrete_angles_0[0]
#print "\t2: ",'%.18e' % discrete_angles_0[1]
#print "\t--- dicrete weights ",energy_0," ---"
#print "\t1: ",'%.18e' % discrete_weights_0[0]
#print "\t2: ",'%.18e' % discrete_weights_0[1]

#discrete_angles_1 = native_data.getMomentPreservingElasticDiscreteAngles(energy_1)
#discrete_weights_1 = native_data.getMomentPreservingElasticWeights(energy_1)
#print "\n\t--- discrete angles ",energy_1," ---"
#print "\t1: ",'%.18e' % discrete_angles_1[0]
#print "\t2: ",'%.18e' % discrete_angles_1[1]
#print "\t--- discrete weights ",energy_1," ---"
#print "\t1: ",'%.18e' % discrete_weights_1[0]
#print "\t2: ",'%.18e' % discrete_weights_1[1]

#print "\n\t--- discrete angles ",energy," ---"
#log_interp = numpy.log(1e-4/1e-5)/numpy.log(1e-3/1e-5)
#lin_interp = (1e-4-1e-5)/(1e-3-1e-5)
#interp = log_interp
#lower_sample = discrete_angles_0[0]

#angle_1 = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*interp
#angle_2 = discrete_angles_0[0] + (discrete_angles_1[1]-discrete_angles_0[0])*interp
#angle_22 = discrete_angles_1[0] + (discrete_angles_1[1]-discrete_angles_0[0])*interp
#angle_3 = discrete_angles_0[1] + (discrete_angles_1[1]-discrete_angles_0[1])*interp

#print "\t1: ",'%.18e' % angle_1
#print "\t2: ",'%.18e' % angle_2
#print "\t3: ",'%.18e' % angle_3

#weight_1_ratio = (discrete_weights_1[0] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\n\tweight_1 cdf range = ','%.16e' % sampling_ratio,' < random number <= ','%.16e' % weight_1_ratio
#weight_2_ratio = (discrete_weights_0[0] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\tweight_2 cdf range = ','%.16e' % weight_1_ratio,' < random number <= ','%.16e' % weight_2_ratio
#print '\tweight_3 cdf range = ','%.16e' % weight_2_ratio,' < random number <= ','%.16e' % 1.0

#random_numbers = [ discrete_weights_1[0]-1e-7, discrete_weights_1[0]+1e-7,\
#                   discrete_weights_0[0]-1e-7, discrete_weights_0[0]+1e-7 ]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "\tsample[",energy,",",random_numbers[0],"] = ", mp_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[1],"] = ", mp_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[2],"] = ", mp_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[3],"] = ", mp_dist.sample( energy )





#print "\n--- lin-lin-lin Tests ---"

####
####  Evaluate Distribution
####
#print '--- Evaluate At Cutoff ---'
#hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, "LinLinLin", True, 1e-7)
#cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, "LinLinLin", True, 1e-7)

#energies = [1e-4, 1e-3, 1e5 ]

#for energy in energies:
#    ratio = hybrid_dist.getSamplingRatio( energy )
#    unorm_eval = cutoff_dist.evaluate(energy, 0.9)
#    evaluation = unorm_eval*ratio

#    print '\nenergy =', energy
#    print '\tsampling ratio =','%.16e' % ratio
#    print '\tunorm eval =','%.16e' % unorm_eval
#    print '\teval   =','%.16e' % evaluation
#    print '\tresult =','%.16e' % hybrid_dist.evaluate( energy, 0.9 )


####
####  Evaluate PDF
####
#print '\n--- Evaluate PDF At Cutoff ---'

#energies = [1e-4, 1e-3, 1e5 ]

#for energy in energies:
#    ratio = hybrid_dist.getSamplingRatio( energy )
#    unorm_pdf = cutoff_dist.evaluatePDF(energy, 0.9)
#    pdf = unorm_pdf*ratio

#    print '\nenergy =', energy
#    print '\tsampling ratio =','%.16e' % ratio
#    print '\tunorm pdf =','%.16e' % unorm_pdf
#    print '\tpdf   =','%.16e' % pdf
#    print '\tresult =','%.16e' % hybrid_dist.evaluatePDF( energy, 0.9 )

####
####  Evaluate CDF
####
#print '\n--- Evaluate CDF At Cutoff ---'

#energies = [1e-4, 1e-3, 1e5 ]

#for energy in energies:
#    ratio = hybrid_dist.getSamplingRatio( energy )
#    unorm_cdf = cutoff_dist.evaluateCDF(energy, 0.9)
#    cdf = unorm_cdf*ratio

#    print '\nenergy =', energy
#    print '\tsampling ratio =','%.16e' % ratio
#    print '\tunorm cdf =','%.16e' % unorm_cdf
#    print '\tcdf   =','%.16e' % cdf
#    print '\tresult =','%.16e' % hybrid_dist.evaluateCDF( energy, 0.9 )


####
####  Cutoff Distribution
####
#print '\n--- Sampling ---'
#energy = 1e-3
#partial_cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, "LinLinLin", True, 1e-14)
#cutoff_cdf = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy)

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy:
#        print 'energy = ', energy_grid[i]
#        print '\tcutoff_cdf = ','%.16e' % cutoff_cdf
#        print '\tcutoff_cs = ','%.16e' % cutoff_cross_sections[i]
#        print '\tmoment_preserving_cs = ','%.16e' % moment_cross_sections[i-moment_index]
#        cross_section_ratio = cutoff_cross_sections[i]*cutoff_cdf/moment_cross_sections[i-moment_index]
#        print '\tcross section ratio = ','%.16e' % cross_section_ratio
#        sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#        print '\tsampling_ratio = ','%.16e' % sampling_ratio

#print "\n\t--- continuous sampling ---"
#cdf = partial_cutoff_dist.evaluateCDF( energy, 0.1 )
#print '\tcdf at 0.1 for cutoff = ','%.16e' % cdf
#cdf_scaled = cdf*sampling_ratio
#print '\tcdf at 0.1 for hybrid = ','%.16e' % cdf_scaled
#print '\tcdf at 0.9 for cutoff = ','%.16e' % cutoff_cdf
#print '\tcdf at 0.9 for hybrid = ','%.16e' % sampling_ratio

#discrete_angles = native_data.getMomentPreservingElasticDiscreteAngles(energy)
#discrete_weights = native_data.getMomentPreservingElasticWeights(energy)

#print "\n\t--- discrete sampling ---"
#print '\tdiscrete_weights[0] = ','%.16e' % discrete_weights[0]
#print '\tdiscrete_angles[0] = ','%.16e' % discrete_angles[0]
#weight_0_ratio = (discrete_weights[0] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\tweight_0 cdf range = ','%.16e' % sampling_ratio,' < random number <= ','%.16e' % weight_0_ratio
#print '\tdiscrete_weights[1] = ','%.16e' % discrete_weights[1]
#print '\tdiscrete_angles[1] = ','%.16e' % discrete_angles[1]
#weight_1_ratio = (discrete_weights[0] + discrete_weights[1] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\tweight_1 cdf range = ','%.16e' % weight_0_ratio,' < random number <= ','%.16e' % weight_1_ratio


#energy = 1e-4
#mp_dist = Collision.createMomentPreservingElasticDistribution(native_data, 0.9, "LinLinLin", True, 1e-15)
#cutoff_cdf = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy )

#for i in range(0, cutoff_cs.size ):
#    if energy_grid[i] == energy:
#        print '\nenergy = ', energy_grid[i]
#        print '\tcutoff_cdf = ','%.16e' % cutoff_cdf
#        print '\tcutoff_cs = ','%.16e' % cutoff_cross_sections[i]
#        print '\tmoment_preserving_cs = ','%.16e' % moment_cross_sections[i-moment_index]
#        cross_section_ratio = cutoff_cross_sections[i]*cutoff_cdf/moment_cross_sections[i-moment_index]
#        print '\tcross section ratio = ','%.16e' % cross_section_ratio
#        cross_section_ratio = hybrid_dist.getCrossSectionRatio( energy )
#        print '\tcross section ratio = ','%.16e' % cross_section_ratio
#        sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#        print '\tsampling_ratio = ','%.16e' % sampling_ratio


#print "\n\t--- continuous sampling ---"
#cdf = partial_cutoff_dist.evaluateCDF( energy, 0.1 )
#print '\tcdf at 0.1 for partial cutoff = ','%.16e' % cdf
#cdf_scaled = cdf*sampling_ratio
#print '\tcdf at 0.1 for hybrid = ','%.16e' % cdf_scaled
#print '\tcdf at 0.9 for cutoff = ','%.16e' % cutoff_cdf
#print '\tcdf at 0.9 for hybrid = ','%.16e' % sampling_ratio

#random_numbers = [cdf,cdf_scaled]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "\tsample[",energy,",",random_numbers[0],"] = ", partial_cutoff_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[1],"] = ",hybrid_dist.sample( energy )

#print "\n\t--- discrete sampling ---"

#index = 0
#for i in range(0, discrete_energy_grid.size ):
#    if discrete_energy_grid[i] <= energy:
#        index = i
#energy_0 = discrete_energy_grid[index]
#energy_1 = discrete_energy_grid[index+1]

#discrete_angles_0 = native_data.getMomentPreservingElasticDiscreteAngles(energy_0)
#discrete_weights_0 = native_data.getMomentPreservingElasticWeights(energy_0)
#print "\n\t--- discrete angles ",energy_0," ---"
#print "\t1: ",'%.18e' % discrete_angles_0[0]
#print "\t2: ",'%.18e' % discrete_angles_0[1]
#print "\t--- dicrete weights ",energy_0," ---"
#print "\t1: ",'%.18e' % discrete_weights_0[0]
#print "\t2: ",'%.18e' % discrete_weights_0[1]

#discrete_angles_1 = native_data.getMomentPreservingElasticDiscreteAngles(energy_1)
#discrete_weights_1 = native_data.getMomentPreservingElasticWeights(energy_1)
#print "\n\t--- discrete angles ",energy_1," ---"
#print "\t1: ",'%.18e' % discrete_angles_1[0]
#print "\t2: ",'%.18e' % discrete_angles_1[1]
#print "\t--- discrete weights ",energy_1," ---"
#print "\t1: ",'%.18e' % discrete_weights_1[0]
#print "\t2: ",'%.18e' % discrete_weights_1[1]

#print "\n\t--- discrete angles ",energy," ---"
#log_interp = numpy.log(1e-4/1e-5)/numpy.log(1e-3/1e-5)
#lin_interp = (1e-4-1e-5)/(1e-3-1e-5)
#interp = lin_interp
#lower_sample = discrete_angles_0[0]

#angle_1 = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*interp
#angle_2 = discrete_angles_0[0] + (discrete_angles_1[1]-discrete_angles_0[0])*interp
#angle_22 = discrete_angles_1[0] + (discrete_angles_1[1]-discrete_angles_0[0])*interp
#angle_3 = discrete_angles_0[1] + (discrete_angles_1[1]-discrete_angles_0[1])*interp

#print "\t1: ",'%.18e' % angle_1
#print "\t2: ",'%.18e' % angle_2
#print "\t3: ",'%.18e' % angle_3

#weight_1_ratio = (discrete_weights_1[0] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\n\tweight_1 cdf range = ','%.16e' % sampling_ratio,' < random number <= ','%.16e' % weight_1_ratio
#weight_2_ratio = (discrete_weights_0[0] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\tweight_2 cdf range = ','%.16e' % weight_1_ratio,' < random number <= ','%.16e' % weight_2_ratio
#print '\tweight_3 cdf range = ','%.16e' % weight_2_ratio,' < random number <= ','%.16e' % 1.0

#random_numbers = [ discrete_weights_1[0]-1e-7, discrete_weights_1[0]+1e-7,\
#                   discrete_weights_0[0]-1e-7, discrete_weights_0[0]+1e-7 ]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "\tsample[",energy,",",random_numbers[0],"] = ", mp_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[1],"] = ", mp_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[2],"] = ", mp_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[3],"] = ", mp_dist.sample( energy )

#### -------------------------------------------------------------------------- ##
####  Adjoint Elastic Unit Test Check
#### -------------------------------------------------------------------------- ##
##print "\n--- Adjoint Elastic Unit Test Check ---"
##data_list = cs_list.get( 'H-Native' )
##native_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
##native_data = Native.AdjointElectronPhotonRelaxationDataContainer( native_file_name )
##energy_grid = native_data.getAdjointElectronEnergyGrid()

##tot_elastic_cs = native_data.getAdjointTotalElasticCrossSection()
##cutoff_cs = native_data.getAdjointCutoffElasticCrossSection()
##screen_rutherford_cs = native_data.getAdjointScreenedRutherfordElasticCrossSection()
##screen_rutherford_index = native_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
##moment_cs = native_data.getAdjointMomentPreservingCrossSection()
##moment_index = native_data.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex()


##discrete_energy_grid = native_data.getAdjointElasticAngularEnergyGrid()

#####
#####  Cutoff Distribution
#####
##energy = 1e-3
##partial_cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, True, True, 1e-15)
##hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, True, True, 1e-15)
##cutoff_cdf = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy )

##for i in range(0, cutoff_cs.size ):
##    if energy_grid[i] <= energy:
##        index = i

##energy_0 = energy_grid[index]
##cutoff_cdf_0 = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy_0 )
##print '\nenergy_0 = ', energy_0
##print '\tcutoff_cdf_0 = ','%.16e' % cutoff_cdf_0
##cutoff_cs_0 = cutoff_cross_sections[index]
##print '\tcutoff_cs_0 = ','%.16e' % cutoff_cs_0
##moment_cs_0 = moment_cross_sections[index-moment_index]
##print '\tmoment_preserving_cs_0 = ','%.16e' % moment_cs_0
##cross_section_ratio_0 = cutoff_cross_sections[index]*cutoff_cdf_0/moment_cross_sections[index-moment_index]
##print '\tcross section ratio_0 = ','%.16e' % cross_section_ratio_0
##cross_section_ratio_0 = hybrid_dist.getCrossSectionRatio( energy_0 )
##print '\tcross section ratio_0 = ','%.16e' % cross_section_ratio_0
##sampling_ratio_0 = cross_section_ratio_0/(1.0+cross_section_ratio_0)
##print '\tsampling_ratio_0 = ','%.16e' % sampling_ratio_0

##energy_1 = energy_grid[index+1]
##cutoff_cdf_1 = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy_1 )
##print '\nenergy_1 = ', energy_1
##print '\tcutoff_cdf_1 = ','%.16e' % cutoff_cdf_1
##cutoff_cs_1 = cutoff_cross_sections[index+1]
##print '\tcutoff_cs_1 = ','%.16e' % cutoff_cs_0
##moment_cs_1 = moment_cross_sections[index+1-moment_index]
##print '\tmoment_preserving_cs_1 = ','%.16e' % moment_cs_1
##cross_section_ratio_1 = cutoff_cross_sections[index+1]*cutoff_cdf_1/moment_cross_sections[index+1-moment_index]
##print '\tcross section ratio_1 = ','%.16e' % cross_section_ratio_1
##cross_section_ratio_1 = hybrid_dist.getCrossSectionRatio( energy_1 )
##print '\tcross section ratio_1 = ','%.16e' % cross_section_ratio_1
##sampling_ratio_1 = cross_section_ratio_1/(1.0+cross_section_ratio_1)
##print '\tsampling_ratio_1 = ','%.16e' % sampling_ratio_1

##lin_interp = (energy - energy_0)/(energy_1 - energy_0)
##log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
##cutoff_cdf = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy )
##print '\nenergy = ', energy
##print '\tcutoff_cdf = ','%.16e' % cutoff_cdf
##cutoff_cs_e = cutoff_cs_0 + (cutoff_cs_1 - cutoff_cs_0)*log_interp
##print '\tcutoff_cs = ','%.16e' % cutoff_cs_e
##moment_cs_e = moment_cs_0 + (moment_cs_1 - moment_cs_0)*log_interp
##print '\tmoment_preserving_cs = ','%.16e' % moment_cs_e
##cross_section_ratio = cutoff_cs_e*cutoff_cdf/moment_cs_e
##print '\tcross section ratio = ','%.16e' % cross_section_ratio
##cross_section_ratio = cross_section_ratio_0 + (cross_section_ratio_1 - cross_section_ratio_0)*lin_interp
##print '\tcross section ratio = ','%.16e' % cross_section_ratio
##cross_section_ratio = hybrid_dist.getCrossSectionRatio( energy )
##print '\tcross section ratio = ','%.16e' % cross_section_ratio
##sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
##print '\tsampling_ratio = ','%.16e' % sampling_ratio


##print "\n\t--- continuous sampling ---"
##cdf_0 = partial_cutoff_dist.evaluateCDF( energy_0, 0.1 )
##cdf = partial_cutoff_dist.evaluateCDF( energy, 0.1 )
##cdf_1 = partial_cutoff_dist.evaluateCDF( energy_1, 0.1 )
##print '\tcdf at 0.1 for partial cutoff = ','%.16e' % cdf
##cdf_0_scaled = cdf_0*sampling_ratio_0
##cdf_scaled = cdf*sampling_ratio
##cdf_1_scaled = cdf_1*sampling_ratio_1
##print '\tcdf at 0.1 for hybrid = ','%.16e' % cdf_scaled
##print '\tcdf at 0.9 for cutoff = ','%.16e' % cutoff_cdf
##print '\tcdf at 0.9 for hybrid = ','%.16e' % sampling_ratio

##random_numbers = [cdf,cdf_scaled]
##Prng.RandomNumberGenerator.setFakeStream(random_numbers)
##print "\tsample[",energy,",",random_numbers[0],"] = ", partial_cutoff_dist.sample( energy )
##print "\tsample[",energy,",",random_numbers[1],"] = ",hybrid_dist.sample( energy )

##print "\n\t--- discrete sampling ---"

##index = 0
##for i in range(0, discrete_energy_grid.size ):
##    if discrete_energy_grid[i] == energy:
##        index = i


##discrete_angles = native_data.getAdjointMomentPreservingElasticDiscreteAngles(energy)
##discrete_weights = native_data.getAdjointMomentPreservingElasticWeights(energy)
##print "\n\t--- discrete angles ",discrete_energy_grid[index]," ---"
##print "\t1: ",'%.18e' % discrete_angles[0]
##print "\t2: ",'%.18e' % discrete_angles[1]
##print "\t--- dicrete weights ",discrete_energy_grid[index]," ---"
##print "\t1: ",'%.18e' % discrete_weights[0]
##print "\t2: ",'%.18e' % discrete_weights[1]


##weight_1_ratio = (discrete_weights[0] + cross_section_ratio)/(1.0+cross_section_ratio)
##print '\n\tweight_1 cdf range = ','%.16e' % sampling_ratio,' < random number <= ','%.16e' % weight_1_ratio
##weight_2_ratio = (discrete_weights[0]+discrete_weights[1] + cross_section_ratio)/(1.0+cross_section_ratio)
##print '\tweight_2 cdf range = ','%.16e' % weight_1_ratio,' < random number <= ','%.16e' % weight_2_ratio

##mp_dist = Collision.createMomentPreservingElasticDistribution(native_data, 0.9, False, True, 1e-15)
##random_numbers = [ discrete_weights[0]-1e-7, discrete_weights[0]+1e-7 ]
##Prng.RandomNumberGenerator.setFakeStream(random_numbers)
##print "\tsample[",energy,",",random_numbers[0],"] = ", mp_dist.sample( energy )
##print "\tsample[",energy,",",random_numbers[1],"] = ", mp_dist.sample( energy )
