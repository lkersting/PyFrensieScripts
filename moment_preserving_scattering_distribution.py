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
#  Elastic Data
# -------------------------------------------------------------------------- ##
elements = ['Al-Native']
interps = ["LinLinLog"]
energies = [1e-3, 1e5, 1e-2 ]
angles = []
for z in elements:
    print "\n----------------------------"
    print "-----", z, "Tests -----"
    print "----------------------------"
    data_list = cs_list.get( z )
    file_name = datadir + data_list.get( 'electroatomic_file_path' )
    native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )

    for interp in interps:
        print "\n--- ",interp,"Tests ---"

        dist = Collision.createMomentPreservingElasticDistribution(native_data, 0.9, interp, True, 1e-7 )

        ###
        ###  Moment Preserving Reaction Unit Test Check
        ###
        print "\nevaluate"
        print "------------------------------------------------"

        for energy in energies:
            if energy == 1e-2:
                lower_angles = native_data.getMomentPreservingElasticDiscreteAngles( 8e-3 )
                upper_angles = native_data.getMomentPreservingElasticDiscreteAngles( 1.6e-2 )
                angles = [lower_angles[0], lower_angles[1], upper_angles[0], upper_angles[1] ]
            else:
                angles = native_data.getMomentPreservingElasticDiscreteAngles( energy )

            print "\n---- energy =", energy, "----"
            for angle in angles:
                pdf = dist.evaluate( energy, angle )
                print '\teval(', '%.16e' % angle,') =\t','%.16e' % pdf

        print "\nevaluatePDF"
        print "------------------------------------------------"

        for energy in energies:
            if energy == 1e-2:
                lower_angles = native_data.getMomentPreservingElasticDiscreteAngles( 8e-3 )
                upper_angles = native_data.getMomentPreservingElasticDiscreteAngles( 1.6e-2 )
                angles = [lower_angles[0], lower_angles[1], upper_angles[0], upper_angles[1] ]
            else:
                angles = native_data.getMomentPreservingElasticDiscreteAngles( energy )

            print "\n---- energy =", energy, "----"
            for angle in angles:
                pdf = dist.evaluatePDF( energy, angle )
                print '\tpdf(', '%.16e' % angle,') =\t','%.16e' % pdf

        print "\nevaluateCDF"
        print "------------------------------------------------"

        for energy in energies:
            if energy == 1e-2:
                lower_angles = native_data.getMomentPreservingElasticDiscreteAngles( 8e-3 )
                upper_angles = native_data.getMomentPreservingElasticDiscreteAngles( 1.6e-2 )
                angles = [lower_angles[0], lower_angles[1], upper_angles[0], upper_angles[1] ]
            else:
                angles = native_data.getMomentPreservingElasticDiscreteAngles( energy )

            print "\n---- energy =", energy, "----"
            for angle in angles:
                cdf = dist.evaluateCDF( energy, angle )
                print '\tcdf(', '%.16e' % angle,') =\t','%.16e' % cdf


        print "\nsample"
        print "------------------------------------------------"

        random_numbers = [ 0.5, 0.5 ]
        for energy in energies:
            Prng.RandomNumberGenerator.setFakeStream(random_numbers)
            print "\n---- energy =", energy, "----"
            e_out, mu = dist.sample( energy )
            print '\tmu =\t','%.16e' % mu
            print '\te_out =\t','%.16e' % e_out



