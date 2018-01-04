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
##  Adjoint Electroionization Data
## -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'C-Native' )
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()

###
###  Electroionization Reaction Unit Test Check
###
print "\n----- C -----"
subshells = adjoint_data.getSubshells()
shell = subshells[0]
adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(shell)

print "\nshell = ", shell

energy = 1e-5
print "energy = ", energy
print '\tcs = ','%.16e' % adjoint_ionization_cs[0]

energy = 1.5
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] <= energy:
        index = i

energy_0 = adjoint_energy_grid[index]
cs_0 = adjoint_ionization_cs[index]
energy_1 = adjoint_energy_grid[index+1]
cs_1 = adjoint_ionization_cs[index+1]

print "energy = ", energy
cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
print '\tcs = ','%.16e' % cs

energy = 1e-3
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] <= energy:
        index = i

energy_0 = adjoint_energy_grid[index]
cs_0 = adjoint_ionization_cs[index]
energy_1 = adjoint_energy_grid[index+1]
cs_1 = adjoint_ionization_cs[index+1]

print "energy = ", energy
cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
print '\tcs = ','%.16e' % cs

energy = 20.0
print "energy = ", energy
print '\tcs = ','%.16e' % adjoint_ionization_cs[adjoint_ionization_cs.size -1]


# Last subshell
shell = subshells[len(subshells) -1]
adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(shell)

print "\nshell = ", shell
energy = 1e-5
print "energy = ", energy
print '\tcs = ','%.16e' % adjoint_ionization_cs[0]

energy = 1.5
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] <= energy:
        index = i

energy_0 = adjoint_energy_grid[index]
cs_0 = adjoint_ionization_cs[index]
energy_1 = adjoint_energy_grid[index+1]
cs_1 = adjoint_ionization_cs[index+1]

print "energy = ", energy
cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
print '\tcs = ','%.16e' % cs

energy = 1e-3
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] <= energy:
        index = i

energy_0 = adjoint_energy_grid[index]
cs_0 = adjoint_ionization_cs[index]
energy_1 = adjoint_energy_grid[index+1]
cs_1 = adjoint_ionization_cs[index+1]

print "energy = ", energy
cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
print '\tcs = ','%.16e' % cs

energy = 20.0
print "energy = ", energy
print '\tcs = ','%.16e' % adjoint_ionization_cs[adjoint_ionization_cs.size -1]

###
###  Electroionization Distribution Unit Test Check
###
print "\n----- H -----"
data_list = cs_list.get( 'H-Native' )
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()
adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(1)


subshells = adjoint_data.getSubshells()
shell = subshells[0]
binding_energy = adjoint_data.getSubshellBindingEnergy( shell )

adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(shell)

print "\nshell = ", shell

energy = 1e-5
print "energy = ", energy
print '\tcs = ','%.16e' % adjoint_ionization_cs[0]

energy = 1e-3
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] <= energy:
        index = i

energy_0 = adjoint_energy_grid[index]
cs_0 = adjoint_ionization_cs[index]
energy_1 = adjoint_energy_grid[index+1]
cs_1 = adjoint_ionization_cs[index+1]

print "energy = ", energy
cs = cs_0 + (cs_1 - cs_0)*( energy - energy_0 )/( energy_1 - energy_0 )
print '\tcs = ','%.16e' % cs

energy = 20.0
print "energy = ", energy
print '\tcs = ','%.16e' % adjoint_ionization_cs[adjoint_ionization_cs.size -1]




interpolations = ["LogLogLog", "LinLinLog" ]

for interp in interpolations:
    ionization_dist = Collision.createLogLogLogUnitBaseCorrelatedElectroionizationSubshellDistribution( adjoint_data, shell, binding_energy, 1e-10)
    if interp == "LinLinLog":
      ionization_dist = Collision.createLinLinLogUnitBaseCorrelatedElectroionizationSubshellDistribution( adjoint_data, shell, binding_energy, 1e-10)
    if interp == "LinLinLin":
      ionization_dist = Collision.createLinLinLinUnitBaseCorrelatedElectroionizationSubshellDistribution( adjoint_data, shell, binding_energy, 1e-10)
    print "\n--- ",interp," ---"
    print "\nEvaluate"

    pdf = ionization_dist.evaluate( 9.99e-6, 2.3711E-5 )
    print "\teval 1 = ",'%.16e' % pdf

    pdf = ionization_dist.evaluate( 1e-5, 2.3711E-5 )
    print "\teval 2 = ",'%.16e' % pdf

    pdf = ionization_dist.evaluate( 1.1e-5, 0.2 )
    print "\teval 3 = ",'%.16e' % pdf

    pdf = ionization_dist.evaluate( 20.0, 20.00002722 )
    print "\teval 4 = ",'%.16e' % pdf

    pdf = ionization_dist.evaluate( 20.01, 22.1 )
    print "\teval 5 = ",'%.16e' % pdf

    print "\nEvaluate PDF"

    pdf = ionization_dist.evaluatePDF( 9.99e-6, 2.3711E-5 )
    print "\tpdf 1 = ",'%.16e' % pdf

    pdf = ionization_dist.evaluatePDF( 1e-5, 2.3711E-5 )
    print "\tpdf 2 = ",'%.16e' % pdf

    pdf = ionization_dist.evaluatePDF( 1.1e-5, 0.2 )
    print "\tpdf 3 = ",'%.16e' % pdf

    pdf = ionization_dist.evaluatePDF( 20.0, 20.00002722 )
    print "\tpdf 4 = ",'%.16e' % pdf

    pdf = ionization_dist.evaluatePDF( 20.01, 22.1 )
    print "\tpdf 5 = ",'%.16e' % pdf

    print "\nEvaluate CDF"

    cdf = ionization_dist.evaluateCDF( 9.99e-6, 1.361e-5 )
    print "\tcdf 1 = ",'%.16e' % cdf

    cdf = ionization_dist.evaluateCDF( 1e-5, 0.2 )
    print "\tcdf 2 = ",'%.16e' % cdf

    cdf = ionization_dist.evaluateCDF( 1.1e-5, 0.2 )
    print "\tcdf 3 = ",'%.16e' % cdf

    cdf = ionization_dist.evaluateCDF( 20.0, 20.00002722 )
    print "\tcdf 4 = ",'%.16e' % cdf

    cdf = ionization_dist.evaluateCDF( 20.01, 22.1 )
    print "\tcdf 5 = ",'%.16e' % cdf


    random_numbers = [ 0.0, 0.1 ]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    incoming_energies = [1e-3, 1e-5]

    i = 0
    for incoming_energy in incoming_energies:
        outgoing_energy, scattering_angle = ionization_dist.sample( incoming_energy )
        print "\n\tincoming_energy =", incoming_energy,"\trandom_number =",random_numbers[i]
        print "\noutgoing_energy = ",'%.16e' % outgoing_energy
        print "scattering_angle = ",'%.16e' % scattering_angle
        i = i+1

###
###  Adjoint Electron Collision Handler Unit Test Check
###
data_list = cs_list.get( 'Si-Native' )
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()

print "\n----- Si -----"
subshells = adjoint_data.getSubshells()
shell = subshells[len(subshells)-1]
adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(shell)

print "\nshell = ", shell

energy = 1.55
index = 0
for i in range(0, adjoint_energy_grid.size ):
    if adjoint_energy_grid[i] <= energy:
        index = i

recoil_energy_0 = adjoint_data.getAdjointElectroionizationRecoilEnergy(shell, adjoint_energy_grid[index])
recoil_pdf_0 = adjoint_data.getAdjointElectroionizationRecoilPDF(shell, adjoint_energy_grid[index])
recoil_energy_1 = adjoint_data.getAdjointElectroionizationRecoilEnergy(shell, adjoint_energy_grid[index+1])
recoil_pdf_1 = adjoint_data.getAdjointElectroionizationRecoilPDF(shell, adjoint_energy_grid[index+1])

print "energy_in_0 = ", adjoint_energy_grid[index]
print "energy_in_1 = ", adjoint_energy_grid[index+1]
print "energy_out_0 = ", recoil_energy_0[0]
print "energy_out_1 = ", recoil_energy_1[0]
print 'pdf_0 = ', recoil_pdf_0[0]
print 'pdf_1 = ', recoil_pdf_1[0]

