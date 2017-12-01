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
###  Cutoff Distribution Pyfrensie Unit Test Check
###
interpolations = ["LinLinLin", "LinLinLog", "LogLogLog"]
cutoff_angle_cosines = [0.9, 0.999999]
energies = [1e5, 1e-3, 4e-4]
angles = [0.0, 0.9]

interpolations = ["LinLinLog"]
cutoff_angle_cosines = [0.9]
for interp in interpolations:
    print "\n\n\t-----",interp,"-----"
    for cutoff in cutoff_angle_cosines:

        cutoff_dist = Collision.createCutoffElasticDistribution(native_data, cutoff, interp, True, 1e-14)
        full_cutoff_dist = Collision.createCutoffElasticDistribution( native_data, 1.0, interp, True, 1e-14)
        print "\n\t--- Cutoff Angle Cosine = ",cutoff," ---"
        print "\n\tEvaluate"
        for energy in energies:
            print "Energy = ",energy
            for angle in angles:
                pdf = cutoff_dist.evaluate( energy, angle )
                full_pdf = full_cutoff_dist.evaluate( energy, angle )
                print '\teval[',angle,']      = ','%.16e' % pdf
#                print '\tfull eval[',angle,'] = ','%.16e' % full_pdf


        print "\n\tEvaluate PDF"
        for energy in energies:
            print "Energy = ",energy
            for angle in angles:
                pdf = cutoff_dist.evaluatePDF( energy, angle )
                full_pdf = full_cutoff_dist.evaluatePDF( energy, angle )
                print '\tPDF[',angle,']      = ','%.16e' % pdf
#                print '\tfull PDF[',angle,'] = ','%.16e' % full_pdf

        print "\n\tEvaluate CDF"
        for energy in energies:
            print "Energy = ",energy
            for angle in angles:
                cdf = cutoff_dist.evaluateCDF( energy, angle )
                full_cdf = full_cutoff_dist.evaluateCDF( energy, angle )
                print '\tCDF[',angle,']      = ','%.16e' % cdf
#                print '\tfull CDF[',angle,'] = ','%.16e' % full_cdf


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

# Get angles and PDF
energies = [1e-5, 1e-3, 1.0e+5]
angular_energy_grid = native_data.getElasticAngularEnergyGrid()
print angular_energy_grid
for energy in energies:
    angles = native_data.getCutoffElasticAngles(energy)
    pdf = native_data.getCutoffElasticPDF(energy)

    index_0 = 0
    for i in range(0, angles.size ):
        if angles[i] <= 0.9:
            index_0 = i

    pdf_at_cutoff = pdf[index_0] + (pdf[index_0+1] - pdf[index_0])*(0.9- angles[index_0] )/(angles[index_0+1] - angles[index_0] );

    print "\n\tgetAngularGridAndPDF Test at energy ",energy

    print 'angle[0]   = ','%.16e' %  angles[0],  '\tpdf[0]     = ','%.16e' %  pdf[0]
    print 'angle[i]   = ','%.16e' %  angles[index_0],  '\tpdf[i]     = ','%.16e' %  pdf[index_0]
    print 'angle      = ','%.16e' %  0.9,  '\tpdf        = ','%.16e' %  pdf_at_cutoff
    print 'angle[i+1] = ','%.16e' %  angles[index_0+1],'\tpdf[i+1]   = ','%.16e' %  pdf[index_0+1]
    print 'angle[N]   = ','%.16e' %  angles[angles.size-1],  '\tpdf[N]     = ','%.16e' %  pdf[angles.size-1]




