#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Data.ACE as ACE
import PyFrensie.Data.ENDL as ENDL
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.Utility.Interpolation as Interpolation
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
import numpy
import matplotlib.pyplot as plt

datadir = '/home/software/mcnpdata/'
#datadir = '/home/lkersting/frensie/src/packages/test_files/'

source = Teuchos.FileInputSource( datadir + '/cross_sections.xml' )
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList( xml_obj )

energy = 15.7

# -------------------------------------------------------------------------- ##
#  Native Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'Au-Native' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
energy_grid = native_data.getElectronEnergyGrid()

index = 0
for i in range(0, energy_grid.size ):
    if energy_grid[i] <= energy:
        index = i

energy_0 = energy_grid[index]
energy_1 = energy_grid[index+1]
lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

# -------------------------------------------------------------------------- ##
#  Native Electroionization Data
# -------------------------------------------------------------------------- ##

total_cs = 0
print "\n----- Au Native -----"
print 'energy_0 = ','%.16e' % energy_0,'\tenergy_1 = ','%.16e' % energy_1


subshells = native_data.getSubshells()
native_cs = [None] * len(subshells)
for i in range(0,len(subshells)):
    shell = subshells[i]
    ionization_cs = native_data.getElectroionizationCrossSection(shell)
    ionization_index = native_data.getElectroionizationCrossSectionThresholdEnergyIndex(shell)

    cs_0 = ionization_cs[index-ionization_index]
    cs_1 = ionization_cs[index+1-ionization_index]
    native_cs[i] = cs_0 + (cs_1 - cs_0)*lin_interp
    #native_cs[i] = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)
    print 'shell = ', shell,'\tcs = ','%.16e' % native_cs[i]

    total_cs += native_cs[i]
print '---------------------------------------------------'
print '\ttotal Native = ','%.16e' % total_cs

# -------------------------------------------------------------------------- ##
#  ACE Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'Au' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
table_name = data_list.get( 'electroatomic_table_name' )
table_start = data_list.get( 'electroatomic_file_start_line' )
ace_file = ACE.ACEFileHandler( file_name, table_name, table_start )
xss_extractor = ACE.XSSEPRDataExtractor( ace_file.getTableNXSArray(), ace_file.getTableJXSArray(), ace_file.getTableXSSArray() )

energy_grid = xss_extractor.extractElectronEnergyGrid()
index = 0
for i in range(0, energy_grid.size ):
    if energy_grid[i] <= energy:
        index = i

energy_0 = energy_grid[index]
energy_1 = energy_grid[index+1]
lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

# -------------------------------------------------------------------------- ##
#  ACE Electroionization Data
# -------------------------------------------------------------------------- ##
raw_subshell_cross_sections = xss_extractor.extractElectroionizationSubshellCrossSections()
subshells = xss_extractor.extractSubshellENDFDesignators()
num_subshells = subshells.size
num_energy_points = energy_grid.size
eion_block = xss_extractor.extractEIONBlock()
eion_loc = xss_extractor.returnEIONLoc()
num_tables = eion_block[0:num_subshells]
table_info = eion_block[num_subshells:2*num_subshells]
table_loc = eion_block[2*num_subshells:3*num_subshells]

ace_cs = [None] * num_subshells
total_cs_ace = 0
print "\n----- Au ACE -----"
print 'energy_0 = ','%.16e' % energy_0,'\tenergy_1 = ','%.16e' % energy_1
for shell in range(0,num_subshells):
    subshell_info = table_info[shell]- eion_loc - 1
    subshell_loc = table_loc[shell]- eion_loc - 1
    ionization_cs = raw_subshell_cross_sections[shell*num_energy_points:(shell+1)*num_energy_points]

    cs_0 = ionization_cs[index]
    cs_1 = ionization_cs[index+1]
    ace_cs[shell] = cs_0 + (cs_1 - cs_0)*lin_interp
    #ace_cs[shell] = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)
    print 'shell = ', subshells[shell],'\tcs = ','%.16e' % ace_cs[shell]

    total_cs_ace += ace_cs[shell]
print '---------------------------------------------------'
print '\ttotal ACE  = ','%.16e' % total_cs_ace

# -------------------------------------------------------------------------- ##
#  ENDL Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'Au-ENDL' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
endl_data = ENDL.ENDLDataContainer( file_name )


# -------------------------------------------------------------------------- ##
#  ENDL Electroionization Data
# -------------------------------------------------------------------------- ##
subshells = endl_data.getSubshells()

total_endl_cs = 0
print "\n----- ENDL Native -----"
endl_cs = [None] * len(subshells)
for i in range(0,len(subshells)):
    shell = subshells[i]
    ionization_cs = endl_data.getElectroionizationCrossSection(shell)
    energy_grid = endl_data.getElectroionizationCrossSectionEnergyGrid(shell)

    index = 0
    for n in range(0, energy_grid.size ):
        if energy_grid[n] <= energy:
            index = n

    energy_0 = energy_grid[index]
    energy_1 = energy_grid[index+1]

    lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
    log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

    cs_0 = ionization_cs[index]
    cs_1 = ionization_cs[index+1]

    #endl_cs[i] = cs_0 + (cs_1 - cs_0)*lin_interp
    endl_cs[i] = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)
    print 'shell = ', shell,'\tcs = ','%.16e' % endl_cs[i]

    total_endl_cs += endl_cs[i]
print '---------------------------------------------------'
print '\ttotal ENDL = ','%.16e' % total_endl_cs

# -------------------------------------------------------------------------- ##
#  Relative Difference in Electroionization Data
# -------------------------------------------------------------------------- ##

print "\n----- Au Relative Difference (Native v ENDL)-----"
for i in range(0,num_subshells):
    rel_diff = (endl_cs[i] - native_cs[i])/endl_cs[i]
    print 'shell = ', subshells[i],'\trel diff = ','%.16e' % rel_diff
print '---------------------------------------------------'
total_rel_diff = (total_endl_cs - total_cs)/total_endl_cs
print '\ttotal rel diff   = ','%.16e' % total_rel_diff

print "\n----- Au Relative Difference (ACE v ENDL)-----"
for i in range(0,num_subshells):
    rel_diff = (endl_cs[i] - ace_cs[i])/endl_cs[i]
    print 'shell = ', subshells[i],'\trel diff = ','%.16e' % rel_diff
print '---------------------------------------------------'
total_rel_diff = (total_endl_cs - total_cs_ace)/total_endl_cs
print '\ttotal rel diff   = ','%.16e' % total_rel_diff


# -------------------------------------------------------------------------- ##
#  Electroionization Data for 1st subshell
# -------------------------------------------------------------------------- ##

# Native
ionization_native_cs = native_data.getElectroionizationCrossSection(subshells[0])
ionization_native_index = native_data.getElectroionizationCrossSectionThresholdEnergyIndex(subshells[0])
native_energy_grid = native_data.getElectronEnergyGrid()
print "\nNative first non zero cross section"
print "\t cross section[",'%.6e' % ionization_native_cs[0],",",'%.6e' % ionization_native_cs[1],"]"
print "\t        energy[",'%.6e' % native_energy_grid[ionization_native_index],",",'%.6e' % native_energy_grid[ionization_native_index+1],"]"


# ACE
ionization_ace_cs = raw_subshell_cross_sections[0:num_energy_points]
ace_energy_grid = xss_extractor.extractElectronEnergyGrid()

index = 0
for n in range(0, ace_energy_grid.size ):
    if ionization_ace_cs[n] == 0.0:
            index = n

print "\nACE first non zero cross section"
print "\t cross section[",'%.6e' % ionization_ace_cs[index],",",'%.6e' % ionization_ace_cs[index+1],"]"
print "\t        energy[",'%.6e' % ace_energy_grid[index],",",'%.6e' % ace_energy_grid[index+1],"]"

# ENDL
ionization_endl_cs = endl_data.getElectroionizationCrossSection(subshells[0])
endl_energy_grid = endl_data.getElectroionizationCrossSectionEnergyGrid(subshells[0])

print "\nENDL first non zero cross section"
print "\t cross section[",'%.6e' % ionization_endl_cs[0],",",'%.6e' % ionization_endl_cs[1],"]"
print "\t        energy[",'%.6e' % endl_energy_grid[0],",",'%.6e' % endl_energy_grid[1],"]"

energy = native_energy_grid[ionization_native_index+1]
energy_0 = endl_energy_grid[0]
energy_1 = endl_energy_grid[1]

lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

cs_0 = ionization_endl_cs[0]
cs_1 = ionization_endl_cs[1]

endl_cs = cs_0 + (cs_1 - cs_0)*lin_interp
#endl_cs = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)
print "energy = ",'%.6e' % energy,"\t cs = ",'%.6e' % endl_cs

energy = ace_energy_grid[index+1]
energy_0 = endl_energy_grid[0]
energy_1 = endl_energy_grid[1]

lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

cs_0 = ionization_endl_cs[0]
cs_1 = ionization_endl_cs[1]

endl_cs = cs_0 + (cs_1 - cs_0)*lin_interp
#endl_cs = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)
print "energy = ",'%.6e' % energy,"\t cs = ",'%.6e' % endl_cs
