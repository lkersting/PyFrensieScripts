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

### -------------------------------------------------------------------------- ##
###  Forward Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'Pb-Native' )
native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )
energy_grid = native_data.getElectronEnergyGrid()

tot_elastic_cs = native_data.getTotalElasticCrossSection()
cutoff_cs = native_data.getCutoffElasticCrossSection()
screen_rutherford_cs = native_data.getScreenedRutherfordElasticCrossSection()
screen_rutherford_index = native_data.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
moment_cs = native_data.getMomentPreservingCrossSection()
moment_index = native_data.getMomentPreservingCrossSectionThresholdEnergyIndex()
discrete_energy_grid = native_data.getElasticAngularEnergyGrid()


interps = ["LogLogLog", "LinLinLin", "LinLinLog"]
energies = [1e-4, 1e-3, 1e5 ]

for interp in interps:
    print "\n----------------------------"
    print "--- ",interp," Pb Tests ---"
    print "----------------------------"

    hybrid_dist = Collision.createLinLinLogExactHybridElasticDistribution(native_data, 0.9, 1e-14)
    cutoff_dist = Collision.createLinLinLogExactCutoffElasticDistribution(native_data, 0.9, 1e-14)
    full_cutoff_dist = Collision.createLinLinLogExactCutoffElasticDistribution(native_data, 1.0, 1e-14)
    if interp == "LinLinLin":
      hybrid_dist = Collision.createLinLinLinExactHybridElasticDistribution(native_data, 0.9, 1e-14)
      cutoff_dist = Collision.createLinLinLinExactCutoffElasticDistribution(native_data, 0.9, 1e-14)
      full_cutoff_dist = Collision.createLinLinLinExactCutoffElasticDistribution(native_data, 1.0, 1e-14)
    elif interp == "LogLogLog":
      hybrid_dist = Collision.createLogLogLogExactHybridElasticDistribution(native_data, 0.9, 1e-14)
      cutoff_dist = Collision.createLogLogLogExactCutoffElasticDistribution(native_data, 0.9, 1e-14)
      full_cutoff_dist = Collision.createLogLogLogExactCutoffElasticDistribution(native_data, 1.0, 1e-14)

    ###
    ###  Get Sampling Ratios
    ###
    print '\n--- Calculate Sampling Ratios ---'

    sampling_ratio = [None]*len(energies)
    for e in range(0, len(energies) ):
        energy = energies[e]
        for i in range(0, cutoff_cs.size ):
            if energy_grid[i] == energy:
                print '\nenergy = ', energy_grid[i]

                cutoff_cdf = cutoff_dist.evaluateCutoffCrossSectionRatio( energy )
                reduced_cs = cutoff_cs[i]*cutoff_cdf
                sampling_ratio[e] = reduced_cs/(moment_cs[i-moment_index]+reduced_cs)

                print '\tcutoff_cdf = ','%.16e' % cutoff_cdf
                print '\tcutoff_cs = ','%.16e' % cutoff_cs[i]
                print '\tmoment_preserving_cs = ','%.16e' % moment_cs[i-moment_index]
                print '\tsampling_ratio = ','%.16e' % sampling_ratio[e]

    ###
    ###  Evaluate Distribution
    ###
    print '\n--- Evaluate At Cutoff ---'

    for i in range(0, len(energies) ):
        energy = energies[i]
        ratio = sampling_ratio[i]
        unorm_eval = cutoff_dist.evaluate(energy, 0.9)
        evaluation = unorm_eval*ratio

        print '\nenergy =', energy
        print '\tsampling ratio =','%.16e' % ratio
        print '\tunorm eval =','%.16e' % unorm_eval
        print '\teval   =','%.16e' % evaluation
        print '\tresult =','%.16e' % hybrid_dist.evaluate( energy, 0.9 )


    ###
    ###  Evaluate PDF
    ###
    print '\n--- Evaluate PDF At Cutoff ---'

    for i in range(0, len(energies) ):
        energy = energies[i]
        ratio = sampling_ratio[i]
        unorm_pdf = cutoff_dist.evaluatePDF(energy, 0.9)
        pdf = unorm_pdf*ratio

        print '\nenergy =', energy
        print '\tsampling ratio =','%.16e' % ratio
        print '\tunorm pdf =','%.16e' % unorm_pdf
        print '\tpdf   =','%.16e' % pdf
        print '\tresult =','%.16e' % hybrid_dist.evaluatePDF( energy, 0.9 )

    ###
    ###  Evaluate CDF
    ###
    print '\n--- Evaluate CDF At Cutoff ---'

    for i in range(0, len(energies) ):
        energy = energies[i]
        ratio = sampling_ratio[i]
        unorm_cdf = cutoff_dist.evaluateCDF(energy, 0.9)
        cdf = unorm_cdf*ratio

        print '\nenergy =', energy
        print '\tsampling ratio =','%.16e' % ratio
        print '\tunorm cdf =','%.16e' % unorm_cdf
        print '\tcdf   =','%.16e' % cdf
        print '\tresult =','%.16e' % hybrid_dist.evaluateCDF( energy, 0.9 )
        print '\tright below =','%.16e' % hybrid_dist.evaluateCDF( energy, 0.9-1e-10 )
        print '\tright above =','%.16e' % hybrid_dist.evaluateCDF( energy, 0.9+1e-10 )


    ###
    ###  Cutoff Distribution
    ###
    print '\n--- Sampling ---'

    energy = 1e-3
    cutoff_cdf = cutoff_dist.evaluateCutoffCrossSectionRatio( energy )
    full_cutoff_cdf = full_cutoff_dist.evaluateCDF( energy, 0.9 )

    for i in range(0, cutoff_cs.size ):
        if energy_grid[i] == energy:
            print '\nenergy = ', energy_grid[i]
            print '\tcutoff_cdf      = ','%.16e' % cutoff_cdf
            print '\tfull_cutoff_cdf = ','%.16e' % full_cutoff_cdf
            print '\tcutoff_cs = ','%.16e' % cutoff_cs[i]
            print '\tmoment_preserving_cs = ','%.16e' % moment_cs[i-moment_index]
            cross_section_ratio = cutoff_cs[i]*cutoff_cdf/moment_cs[i-moment_index]
            print '\tcross section ratio = ','%.16e' % cross_section_ratio
            sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
            print '\tsampling_ratio = ','%.16e' % sampling_ratio
            print "\tsampling_ratio = 0.24463192152870505414"
    print "\n\t--- continuous sampling ---"
    cdf = cutoff_dist.evaluateCDF( energy, 0.1 )
    print '\tcdf at 0.1 for cutoff = ','%.16e' % cdf
    cdf_scaled = cdf*sampling_ratio
    print '\tcdf at 0.1 for hybrid = ','%.16e' % cdf_scaled
    print '\tcdf at 0.9 for cutoff = ','%.16e' % cutoff_cdf
    print '\tcdf at 0.9 for hybrid = ','%.16e' % sampling_ratio

    random_numbers = [cdf,cdf_scaled]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
    print "\tsample[",energy,",",random_numbers[0],"] = ", cutoff_dist.sample( energy )
    print "\tsample[",energy,",",random_numbers[1],"] = ", hybrid_dist.sample( energy )

    discrete_angles = native_data.getMomentPreservingElasticDiscreteAngles(energy)
    discrete_weights = native_data.getMomentPreservingElasticWeights(energy)

    print "\n\t--- discrete sampling ---"
    print '\tdiscrete_weights[0] = ','%.16e' % discrete_weights[0]
    print '\tdiscrete_angles[0] = ','%.16e' % discrete_angles[0]
    weight_0_ratio = (discrete_weights[0] + cross_section_ratio)/(1.0+cross_section_ratio)
    print '\tweight_0 cdf range = ','%.16e' % sampling_ratio,' < random number <= ','%.16e' % weight_0_ratio
    print '\tdiscrete_weights[1] = ','%.16e' % discrete_weights[1]
    print '\tdiscrete_angles[1] = ','%.16e' % discrete_angles[1]
    weight_1_ratio = (discrete_weights[0] + discrete_weights[1] + cross_section_ratio)/(1.0+cross_section_ratio)
    print '\tweight_1 cdf range = ','%.16e' % weight_0_ratio,' < random number <= ','%.16e' % weight_1_ratio


    energy = 1e-4
    mp_dist = Collision.createLinLinLogExactMomentPreservingElasticDistribution(native_data, 0.9, 1e-15)
    if interp == "LinLinLin":
      mp_dist = Collision.createLinLinLinExactMomentPreservingElasticDistribution(native_data, 0.9, 1e-15)
    elif interp == "LogLogLog":
      mp_dist = Collision.createLogLogLogExactMomentPreservingElasticDistribution(native_data, 0.9, 1e-15)
    cutoff_cdf = cutoff_dist.evaluateCutoffCrossSectionRatio( energy )

    for i in range(0, cutoff_cs.size ):
        if energy_grid[i] == energy:
            print '\nenergy = ', energy_grid[i]
            print '\tcutoff_cdf = ','%.16e' % cutoff_cdf
            print '\tcutoff_cs = ','%.16e' % cutoff_cs[i]
            print '\tmoment_preserving_cs = ','%.16e' % moment_cs[i-moment_index]
            cross_section_ratio = cutoff_cs[i]*cutoff_cdf/moment_cs[i-moment_index]
            print '\tcross section ratio = ','%.16e' % cross_section_ratio
            sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
            print '\tsampling_ratio = ','%.16e' % sampling_ratio


    print "\n\t--- continuous sampling ---"
    cdf = cutoff_dist.evaluateCDF( energy, 0.1 )
    print '\tcdf at 0.1 for partial cutoff = ','%.16e' % cdf
    cdf_scaled = cdf*sampling_ratio
    cdf_scaled = hybrid_dist.evaluateCDF( energy, 0.1 )
    print '\tcdf at 0.1 for hybrid = ','%.16e' % cdf_scaled
    print '\tcdf at 0.9 for cutoff = ','%.16e' % cutoff_cdf
    cdf_scaled = hybrid_dist.evaluateCDF( energy, 0.7 )
    print '\tcdf at 0.9 for hybrid = ','%.16e' % cdf_scaled

    random_numbers = [cdf,cdf_scaled]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
    print "\tsample[",energy,",",random_numbers[0],"] = ", cutoff_dist.sample( energy )
    print "\tsample[",energy,",",random_numbers[1],"] = ", hybrid_dist.sample( energy )

    print "\n\t--- discrete sampling ---"

    index = 0
    for i in range(0, discrete_energy_grid.size ):
        if discrete_energy_grid[i] <= energy:
            index = i
    energy_0 = discrete_energy_grid[index]
    energy_1 = discrete_energy_grid[index+1]

    discrete_angles_0 = native_data.getMomentPreservingElasticDiscreteAngles(energy_0)
    discrete_weights_0 = native_data.getMomentPreservingElasticWeights(energy_0)
    print "\n\t--- discrete angles ",energy_0," ---"
    print "\t1: ",'%.18e' % discrete_angles_0[0]
    print "\t2: ",'%.18e' % discrete_angles_0[1]
    print "\t--- dicrete weights ",energy_0," ---"
    print "\t1: ",'%.18e' % discrete_weights_0[0]
    print "\t2: ",'%.18e' % discrete_weights_0[1]

    discrete_angles_1 = native_data.getMomentPreservingElasticDiscreteAngles(energy_1)
    discrete_weights_1 = native_data.getMomentPreservingElasticWeights(energy_1)
    print "\n\t--- discrete angles ",energy_1," ---"
    print "\t1: ",'%.18e' % discrete_angles_1[0]
    print "\t2: ",'%.18e' % discrete_angles_1[1]
    print "\t--- discrete weights ",energy_1," ---"
    print "\t1: ",'%.18e' % discrete_weights_1[0]
    print "\t2: ",'%.18e' % discrete_weights_1[1]

    print "\n\t--- discrete angles ",energy," ---"

    log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
    lin_interp = (energy-energy_0)/(energy_1-energy_0)

    angle_1 = 0.0
    angle_2 = 0.0
    angle_3 = 0.0

    if interp =="LogLogLog":
        print "\tusing log-log interp"
        angle_1 = 1.0-(1.0-discrete_angles_0[0])*pow(((1.0-discrete_angles_1[0])/(1.0-discrete_angles_0[0])),log_interp)
        angle_2 = 1.0-(1.0-discrete_angles_0[0])*pow(((1.0-discrete_angles_1[1])/(1.0-discrete_angles_0[0])),log_interp)
        angle_22 = 1.0-(1.0-discrete_angles_1[0])*pow(((1.0-discrete_angles_1[1])/(1.0-discrete_angles_0[0])),log_interp)
        angle_3 = 1.0-(1.0-discrete_angles_0[1])*pow(((1.0-discrete_angles_1[1])/(1.0-discrete_angles_0[1])),log_interp)
    elif interp == "LinLinLin":
        print "\tusing lin-lin interp"
        angle_1 = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*lin_interp
        angle_2 = discrete_angles_0[0] + (discrete_angles_1[1]-discrete_angles_0[0])*lin_interp
        angle_22 = discrete_angles_1[0] + (discrete_angles_1[1]-discrete_angles_0[0])*lin_interp
        angle_3 = discrete_angles_0[1] + (discrete_angles_1[1]-discrete_angles_0[1])*lin_interp
    else:
        print "\tusing lin-log interp"
        angle_1 = discrete_angles_0[0] + (discrete_angles_1[0]-discrete_angles_0[0])*log_interp
        angle_2 = discrete_angles_0[0] + (discrete_angles_1[1]-discrete_angles_0[0])*log_interp
        angle_22 = discrete_angles_1[0] + (discrete_angles_1[1]-discrete_angles_0[0])*log_interp
        angle_3 = discrete_angles_0[1] + (discrete_angles_1[1]-discrete_angles_0[1])*log_interp

    print "\t1: ",'%.18e' % angle_1
    print "\t2: ",'%.18e' % angle_2
    print "\t3: ",'%.18e' % angle_3

    weight_1_ratio = (discrete_weights_1[0] + cross_section_ratio)/(1.0+cross_section_ratio)
    print '\n\tweight_1 cdf range = ','%.16e' % sampling_ratio,' < random number <= ','%.16e' % weight_1_ratio
    weight_2_ratio = (discrete_weights_0[0] + cross_section_ratio)/(1.0+cross_section_ratio)
    print '\tweight_2 cdf range = ','%.16e' % weight_1_ratio,' < random number <= ','%.16e' % weight_2_ratio
    print '\tweight_3 cdf range = ','%.16e' % weight_2_ratio,' < random number <= ','%.16e' % 1.0

    random_numbers = [ discrete_weights_1[0]-1e-7, discrete_weights_1[0]+1e-7,\
                       discrete_weights_0[0]-1e-7, discrete_weights_0[0]+1e-7 ]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
    print "\tsample[",energy,",",random_numbers[0],"] = ", mp_dist.sample( energy )
    print "\tsample[",energy,",",random_numbers[1],"] = ", mp_dist.sample( energy )
    print "\tsample[",energy,",",random_numbers[2],"] = ", mp_dist.sample( energy )
    print "\tsample[",energy,",",random_numbers[3],"] = ", mp_dist.sample( energy )



## -------------------------------------------------------------------------- ##
##  Adjoint Elastic Unit Test Check
## -------------------------------------------------------------------------- ##
print "\n--- Adjoint Elastic Unit Test Check ---"
data_list = cs_list.get( 'H-Native' )
native_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
native_data = Native.AdjointElectronPhotonRelaxationDataContainer( native_file_name )
energy_grid = native_data.getAdjointElectronEnergyGrid()

tot_elastic_cs = native_data.getAdjointTotalElasticCrossSection()
cutoff_cs = native_data.getAdjointCutoffElasticCrossSection()
screen_rutherford_cs = native_data.getAdjointScreenedRutherfordElasticCrossSection()
screen_rutherford_index = native_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()
moment_cs = native_data.getAdjointMomentPreservingCrossSection()
moment_index = native_data.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex()


discrete_energy_grid = native_data.getAdjointElasticAngularEnergyGrid()

###
###  Cutoff Distribution
###
energy = 1e-3

partial_cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, True, True, 1e-15)
full_cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 1.0, True, True, 1e-15)
hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, True, True, 1e-15)
if interp == "LinLinLin":
  mp_dist = Collision.createLinLinLinExactMomentPreservingElasticDistribution(native_data, 0.9, 1e-15)
elif interp == "LogLogLog":
  mp_dist = Collision.createLogLogLogExactMomentPreservingElasticDistribution(native_data, 0.9, 1e-15)

partial_cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 0.9, True, True, 1e-15)
full_cutoff_dist = Collision.createCutoffElasticDistribution(native_data, 1.0, True, True, 1e-15)
hybrid_dist = Collision.createHybridElasticDistribution(native_data, 0.9, True, True, 1e-15)
cutoff_cdf = full_cutoff_dist.evaluateCDF( energy, 0.9 )

for i in range(0, cutoff_cs.size ):
    if energy_grid[i] <= energy:
        index = i

energy_0 = energy_grid[index]
cutoff_cdf_0 = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy_0 )
print '\nenergy_0 = ', energy_0
print '\tcutoff_cdf_0 = ','%.16e' % cutoff_cdf_0
cutoff_cs_0 = cutoff_cs[index]
print '\tcutoff_cs_0 = ','%.16e' % cutoff_cs_0
moment_cs_0 = moment_cs[index-moment_index]
print '\tmoment_preserving_cs_0 = ','%.16e' % moment_cs_0
cross_section_ratio_0 = cutoff_cs[index]*cutoff_cdf_0/moment_cs[index-moment_index]
print '\tcross section ratio_0 = ','%.16e' % cross_section_ratio_0
#cross_section_ratio_0 = hybrid_dist.getCrossSectionRatio( energy_0 )
#print '\tcross section ratio_0 = ','%.16e' % cross_section_ratio_0
#sampling_ratio_0 = cross_section_ratio_0/(1.0+cross_section_ratio_0)
#print '\tsampling_ratio_0 = ','%.16e' % sampling_ratio_0

energy_1 = energy_grid[index+1]
cutoff_cdf_1 = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy_1 )
print '\nenergy_1 = ', energy_1
print '\tcutoff_cdf_1 = ','%.16e' % cutoff_cdf_1
cutoff_cs_1 = cutoff_cs[index+1]
print '\tcutoff_cs_1 = ','%.16e' % cutoff_cs_0
moment_cs_1 = moment_cs[index+1-moment_index]
print '\tmoment_preserving_cs_1 = ','%.16e' % moment_cs_1
cross_section_ratio_1 = cutoff_cs[index+1]*cutoff_cdf_1/moment_cs[index+1-moment_index]
print '\tcross section ratio_1 = ','%.16e' % cross_section_ratio_1
#cross_section_ratio_1 = hybrid_dist.getCrossSectionRatio( energy_1 )
#print '\tcross section ratio_1 = ','%.16e' % cross_section_ratio_1
#sampling_ratio_1 = cross_section_ratio_1/(1.0+cross_section_ratio_1)
#print '\tsampling_ratio_1 = ','%.16e' % sampling_ratio_1

lin_interp = (energy - energy_0)/(energy_1 - energy_0)
log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
cutoff_cdf = partial_cutoff_dist.evaluateCutoffCrossSectionRatio( energy )
print '\nenergy = ', energy
print '\tcutoff_cdf = ','%.16e' % cutoff_cdf
cutoff_cs_e = cutoff_cs_0 + (cutoff_cs_1 - cutoff_cs_0)*log_interp
print '\tcutoff_cs = ','%.16e' % cutoff_cs_e
moment_cs_e = moment_cs_0 + (moment_cs_1 - moment_cs_0)*log_interp
print '\tmoment_preserving_cs = ','%.16e' % moment_cs_e
cross_section_ratio = cutoff_cs_e*cutoff_cdf/moment_cs_e
print '\tcross section ratio = ','%.16e' % cross_section_ratio
#cross_section_ratio = cross_section_ratio_0 + (cross_section_ratio_1 - cross_section_ratio_0)*lin_interp
#print '\tcross section ratio = ','%.16e' % cross_section_ratio
#cross_section_ratio = hybrid_dist.getCrossSectionRatio( energy )
#print '\tcross section ratio = ','%.16e' % cross_section_ratio
#sampling_ratio = cross_section_ratio/(1.0+cross_section_ratio)
#print '\tsampling_ratio = ','%.16e' % sampling_ratio


print "\n\t--- continuous sampling ---"
cdf_0 = partial_cutoff_dist.evaluateCDF( energy_0, 0.1 )
cdf = partial_cutoff_dist.evaluateCDF( energy, 0.1 )
cdf_1 = partial_cutoff_dist.evaluateCDF( energy_1, 0.1 )
print '\tcdf at 0.1 for partial cutoff = ','%.16e' % cdf
cdf_scaled = cdf*sampling_ratio
print '\tcdf at 0.1 for hybrid = ','%.16e' % cdf_scaled
print '\tcdf at 0.9 for cutoff = ','%.16e' % cutoff_cdf
print '\tcdf at 0.9 for hybrid = ','%.16e' % sampling_ratio

random_numbers = [cdf,cdf_scaled]
Prng.RandomNumberGenerator.setFakeStream(random_numbers)
print "\tsample[",energy,",",random_numbers[0],"] = ", partial_cutoff_dist.sample( energy )
print "\tsample[",energy,",",random_numbers[1],"] = ",hybrid_dist.sample( energy )

print "\n\t--- discrete sampling ---"

index = 0
for i in range(0, discrete_energy_grid.size ):
    if discrete_energy_grid[i] == energy:
        index = i


discrete_angles = native_data.getAdjointMomentPreservingElasticDiscreteAngles(energy)
discrete_weights = native_data.getAdjointMomentPreservingElasticWeights(energy)
print "\n\t--- discrete angles ",discrete_energy_grid[index]," ---"
print "\t1: ",'%.18e' % discrete_angles[0]
print "\t2: ",'%.18e' % discrete_angles[1]
#print "\t--- dicrete weights ",discrete_energy_grid[index]," ---"
#print "\t1: ",'%.18e' % discrete_weights[0]
#print "\t2: ",'%.18e' % discrete_weights[1]


#weight_1_ratio = (discrete_weights[0] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\n\tweight_1 cdf range = ','%.16e' % sampling_ratio,' < random number <= ','%.16e' % weight_1_ratio
#weight_2_ratio = (discrete_weights[0]+discrete_weights[1] + cross_section_ratio)/(1.0+cross_section_ratio)
#print '\tweight_2 cdf range = ','%.16e' % weight_1_ratio,' < random number <= ','%.16e' % weight_2_ratio

#mp_dist = Collision.createMomentPreservingElasticDistribution(native_data, 0.9, False, True, 1e-15)
#random_numbers = [ discrete_weights[0]-1e-7, discrete_weights[0]+1e-7 ]
#Prng.RandomNumberGenerator.setFakeStream(random_numbers)
#print "\tsample[",energy,",",random_numbers[0],"] = ", mp_dist.sample( energy )
#print "\tsample[",energy,",",random_numbers[1],"] = ", mp_dist.sample( energy )

