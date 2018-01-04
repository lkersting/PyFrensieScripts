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
#  Electroionization Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'C-Native' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
energy_grid = native_data.getElectronEnergyGrid()

###
###  Electroionization Reaction Unit Test Check
###
print "\n----- C -----"
subshells = native_data.getSubshells()
shell = subshells[0]
binding_energy = native_data.getSubshellBindingEnergy( shell )

ionization_cs = native_data.getElectroionizationCrossSection(shell)
threshold = native_data.getElectroionizationCrossSectionThresholdEnergyIndex(shell)
ionization_dist = Collision.createLinLinLogUnitBaseCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-15)

print "\nshell = ", shell

energies = [1.70425200079801E-03, 1.70425200079802E-03,1.98284583249127E-03,2.00191878322064E-03]
e_outs = [8.52126000399011E-04, 8.52126000399011E-04, 8.52126000399011E-04, 8.52126000399011E-04]

for i in range(0,len(energies)):
    energy = energies[i]
    e_out = e_outs[i]
    print "\nenergy = ", energy, "\tknock-on energy =",e_out

    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    energy_0 = energy_grid[index]
    cs_0 = ionization_cs[index-threshold]
    energy_1 = energy_grid[index+1]
    cs_1 = ionization_cs[index+1-threshold]

    log_interp = numpy.log( energy/energy_0 )/numpy.log( energy_1/energy_0 )
    cs =  (cs_0)*pow((cs_1/cs_0),log_interp)
    print '\tcs = ','%.16e' % cs
    pdf = ionization_dist.evaluatePDF( energy, e_out )
    print '\tpdf = ','%.16e' % pdf
    dcs = cs*pdf
    print '\tdcs = ','%.16e' % dcs



shell = subshells[len(subshells)-1]
binding_energy = native_data.getSubshellBindingEnergy( shell )

ionization_cs = native_data.getElectroionizationCrossSection(shell)
threshold = native_data.getElectroionizationCrossSectionThresholdEnergyIndex(shell)
ionization_dist = Collision.createLinLinLogUnitBaseCorrelatedElectroionizationSubshellDistribution( native_data, shell, binding_energy, 1e-12)

print "\nshell = ", shell

energies = [0.0025118800000459599528, 0.0025118800000459773,0.002511885,0.0025118897153524992472,0.0025118908794333669708]
e_outs = [0.0012514500000459765489, 0.0012514500000459765489, 0.0012514500000459765489, 0.0012514500000459765489, 0.0012514500000459765489]

for i in range(0,len(energies)):
    energy = energies[i]
    e_out = e_outs[i]
    print "\nenergy = ", energy, "\tknock-on energy =",e_out

    index = 0
    for i in range(0, energy_grid.size ):
        if energy_grid[i] <= energy:
            index = i

    energy_0 = energy_grid[index]
    cs_0 = ionization_cs[index-threshold]
    energy_1 = energy_grid[index+1]
    cs_1 = ionization_cs[index+1-threshold]

    log_interp = numpy.log( energy/energy_0 )/numpy.log( energy_1/energy_0 )
    cs =  (cs_0)*pow((cs_1/cs_0),log_interp)
    print '\tcs = ','%.16e' % cs
    pdf = ionization_dist.evaluatePDF( energy, e_out )
    print '\tpdf = ','%.16e' % pdf
    dcs = cs*pdf
    print '\tdcs = ','%.16e' % dcs