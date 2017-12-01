#! /usr/bin/env python
import PyFrensie.Data.ACE as ACE
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
file_names = [ '/home/lkersting/frensie/src/packages/test_files/ace/test_h_epr14_ace_file.txt']
table_names = ['1000.14p']

for j in range(0,len(table_names) ):
    print "\n----------------------------"
    print "-----", table_names[j], "Tests -----"
    print "----------------------------"
    ace_file = ACE.ACEFileHandler( file_names[j], table_names[j], 1 )
    ace_data = ACE.XSSEPRDataExtractor( ace_file.getTableNXSArray(), ace_file.getTableJXSArray(), ace_file.getTableXSSArray() )

    cutoff_dist = Collision.createCutoffElasticDistribution( ace_data )

    ###
    ###  Cutoff Distribution Pyfrensie Unit Test Check
    ###
    energies = [1e5, 1e-3, 4e-4]
    angles = [0.0, 0.9]

    print "\n\tEvaluate"
    for energy in energies:
        print "Energy = ",energy
        for angle in angles:
            pdf = cutoff_dist.evaluate( energy, angle )
            print '\teval[',angle,']      = ','%.16e' % pdf


    # print "\n\tEvaluate PDF"
    # for energy in energies:
    #     print "Energy = ",energy
    #     for angle in angles:
    #         pdf = cutoff_dist.evaluatePDF( energy, angle )
    #         print '\tPDF[',angle,']      = ','%.16e' % pdf

    print "\n\tEvaluate CDF"
    for energy in energies:
        print "Energy = ",energy
        for angle in angles:
            cdf = cutoff_dist.evaluateCDF( energy, angle )
            print '\tCDF[',angle,']      = ','%.16e' % cdf


    energy = 1e-3
    random_numbers = [ 0.5 ]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    print "\n\t--- Sample ---"
    print "Energy = ",energy
    energy,angle = cutoff_dist.sample( energy )
    print 'angle[0.5]         = ','%.16e' % angle


    energy = 1e5
    random_numbers = [ 0.5]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    print "Energy = ",energy
    energy,angle = cutoff_dist.sample( energy )
    print 'angle[0.5]         = ','%.16e' % angle
