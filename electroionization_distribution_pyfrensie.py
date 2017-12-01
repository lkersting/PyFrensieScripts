#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.Utility.Interpolation as Interpolation
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
import numpy

Utility.initFrensiePrng()

#datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/frensie/src/packages/test_files/'

source = Teuchos.FileInputSource( datadir + '/cross_sections.xml' )
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList( xml_obj )

# -------------------------------------------------------------------------- ##
#  Electroionization Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'H-Native' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
energy_grid = native_data.getElectronEnergyGrid()

###
###  Electroionization Reaction Unit Test Check
###
print "\n----- H -----"
subshells = native_data.getSubshells()
shell = subshells[0]
binding_energy = native_data.getSubshellBindingEnergy( shell )

ionization_dist = Collision.createLinLinLogCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-7)

print "\nshell = ", shell

energies = [1e5, 1e5, 1e5, 1.0e-3, 1.0e-3, 1.0e-3]
e_outs = [1e-5, 1.0, 10.0, 1e-5, 1e-4, 5e-4]

print "\n--- evaluate ---\n";

for i in range(0,len(energies)):
    energy = energies[i]
    e_out = e_outs[i]
    print "energy = ", energy, "\tknock-on energy =",e_out

    pdf = ionization_dist.evaluate( energy, e_out )
    print '\tevaluate = ','%.16e' % pdf

print "\n--- evaluate PDF ---\n";

for i in range(0,len(energies)):
    energy = energies[i]
    e_out = e_outs[i]
    print "energy = ", energy, "\tknock-on energy =",e_out

    pdf = ionization_dist.evaluatePDF( energy, e_out )
    print '\tpdf = ','%.16e' % pdf