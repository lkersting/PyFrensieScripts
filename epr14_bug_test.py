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

# -------------------------------------------------------------------------- ##
#  Electroionization Data
# -------------------------------------------------------------------------- ##
file_names = ['/home/lkersting/frensie/src/packages/test_files/ace/test_au_epr14_ace_file.txt']
table_names = ['79000.14p']

for j in range(0,len(table_names) ):
    print "\n----------------------------"
    print "-----", table_names[j], "Tests -----"
    print "----------------------------"
    ace_file = ACE.ACEFileHandler( file_names[j], table_names[j], 1 )
    ace_data = ACE.XSSEPRDataExtractor( ace_file.getTableNXSArray(), ace_file.getTableJXSArray(), ace_file.getTableXSSArray() )

    subshells = ace_data.extractSubshellENDFDesignators()
    print subshells
    shell_binding_energies = ace_data.extractSubshellBindingEnergies()
    shells = [subshells[0], subshells[len(subshells)-1] ]
    print shells
    binding_energies = [shell_binding_energies[0], shell_binding_energies[len(shell_binding_energies)-1] ]

    shells = [subshells[0]]
    binding_energies = [shell_binding_energies[0]]

    for i in range(0, len(shells) ):
        shell = int(shells[i])
        binding_energy = binding_energies[i]
        print "\nshell = ", shell, "\tbinding energy =", binding_energy;
        print "----------------------------\n"

        ionization_dist = Collision.createElectroionizationSubshellDistribution( ace_data, shell, 1e-7 )

        print "\n--- evaluateCDF ---\n";

        energies = [binding_energy + 1e-8, binding_energy + 3e-8, 9.12175e-2, 1.0, 1.0, 1e5]
        e_outs = [1e-8, 1.0001e-8, 4.275e-4, 1.33136131511529e-1, 9.7163e-2, 1.75297e2]

        energies = [1e5, 1e5, 1e5, 1e-3, 1e-3, 1e-3]
        e_outs = [1e-5, 1.0, 10.0, 1e-5, 1e-4, 5e-4]

        for i in range(0,len(energies)):
            cdf = ionization_dist.evaluateCDF( energies[i], e_outs[i] )
            print "\tevaluateCDF[",'%.6e' % energies[i],",",'%.10e' % e_outs[i],"] =",'%.16e' % cdf


        print "\n--- sample ---\n";

        energies = [1.0, 1.5]
        random_numbers = [0.5, 0.5]

        energies = [1e-3, 1e5]
        random_numbers = [0.0, 1.0-1e-15]
        e_out = 0.0
        knock_on_energy = 0.0
        scattering_angle_cosine = 0.0
        knock_on_angle_cosine = 0.0

        for i in range(0,len(energies)):

            Prng.RandomNumberGenerator.setFakeStream(random_numbers)

            e_out, knock_on_energy, scattering_angle_cosine, knock_on_angle_cosine = ionization_dist.samplePrimaryAndSecondary( energies[i] )
            print "\tenergies[i] =",energies[i],"\trandom_number =",random_numbers[0]
            print "knock-on mu   =",'%.16e' % knock_on_angle_cosine, "\tknock-on energy =",'%.16e' % knock_on_energy
            print "scattering mu =",'%.16e' % scattering_angle_cosine, "\toutgoing energy =",'%.16e' % e_out

            e_out, knock_on_energy, scattering_angle_cosine, knock_on_angle_cosine = ionization_dist.samplePrimaryAndSecondary( energies[i] )
            print "\tenergies[i] =",energies[i],"\trandom_number =",random_numbers[0]
            print "knock-on mu   =",'%.16e' % knock_on_angle_cosine, "\tknock-on energy =",'%.16e' % knock_on_energy
            print "scattering mu =",'%.16e' % scattering_angle_cosine, "\toutgoing energy =",'%.16e' % e_out

