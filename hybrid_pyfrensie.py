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

### -------------------------------------------------------------------------- ##
###  Forward Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )

###
###  Hybrid Distribution PyFrensie Unit Test Check
###
interpolations = ["LinLinLin", "LinLinLog", "LogLogLog"]
energies = [1e5, 1e-3, 4e-4]
angles = [0.0, 0.9]

energies = [1e5,1e-3]
interpolations = ["LinLinLog"]

for interp in interpolations:
    print "\n\n\t-----",interp,"-----"

    hybrid_dist = Collision.createLinLinLogCorrelatedHybridElasticDistribution(native_data, 0.9, 1e-15)

    for energy in energies:
        print "Energy = ",energy

        print "\tEvaluate"
        for angle in angles:
            pdf = hybrid_dist.evaluate( energy, angle )
            print '\teval[',angle,'] = ','%.16e' % pdf


    for energy in energies:
        print "Energy = ",energy
        print "\tEvaluate PDF"
        for angle in angles:
            pdf = hybrid_dist.evaluatePDF( energy, angle )
            print '\tPDF[',angle,'] = ','%.16e' % pdf


    for energy in energies:
        print "Energy = ",energy
        print "\tEvaluate CDF"
        for angle in angles:
            cdf = hybrid_dist.evaluateCDF( energy, angle )
            print '\tCDF[',angle,'] = ','%.16e' % cdf


    energy = 1e-3
    random_numbers = [ 0.5 ]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    print "\n--- Sample ---"
    print "Energy = ",energy
    energy,angle = hybrid_dist.sample( energy )
    print 'angle[0.5]         = ','%.16e' % angle


    energy = 1e5
    random_numbers = [ 0.5]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    print "\nEnergy = ",energy
    energy,angle = hybrid_dist.sample( energy )
    print 'angle[0.5]         = ','%.16e' % angle

for energy in energies:
    print "\nEnergy = ",energy
    angles = native_data.getMomentPreservingElasticDiscreteAngles( energy )
    weights = native_data.getMomentPreservingElasticWeights( energy )
    for i in range(0, len(angles) ):
        print "\tangle[",i,"]  =",'%.16e' % angles[i]
        print "\tweight[",i,"] =",'%.16e' % weights[i]

