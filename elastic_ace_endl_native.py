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
native_energy_grid = native_data.getElectronEnergyGrid()

analog_dist = Collision.createAnalogElasticDistribution(native_data, True, True, 1e-15)
#lin_analog_dist = Collision.createAnalogElasticDistribution(native_data, False, True, 1e-15)

native_index = 0
for i in range(0, native_energy_grid.size ):
    if native_energy_grid[i] <= energy:
        native_index = i

energy_0 = native_energy_grid[native_index]
energy_1 = native_energy_grid[native_index+1]
native_lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
native_log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

# -------------------------------------------------------------------------- ##
#  Native Cutoff Elastic Cross Section
# -------------------------------------------------------------------------- ##
elastic_cross_sections = native_data.getCutoffElasticCrossSection()

cs_0 = elastic_cross_sections[native_index]
cs_1 = elastic_cross_sections[native_index+1]
native_elastic_cs = cs_0 + (cs_1 - cs_0)*native_lin_interp
#native_elastic_cs = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*native_log_interp)

# -------------------------------------------------------------------------- ##
#  Native Analog Elastic Cross Section
# -------------------------------------------------------------------------- ##
native_analog_cs = native_elastic_cs/analog_dist.evaluateCDFAtCutoff( energy )
#native_elastic_cs = native_elastic_cs/lin_analog_dist.evaluateCDFAtCutoff( energy )

# -------------------------------------------------------------------------- ##
#  ACE Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'Au' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
table_name = data_list.get( 'electroatomic_table_name' )
table_start = data_list.get( 'electroatomic_file_start_line' )
ace_file = ACE.ACEFileHandler( file_name, table_name, table_start )
xss_extractor = ACE.XSSEPRDataExtractor( ace_file.getTableNXSArray(), ace_file.getTableJXSArray(), ace_file.getTableXSSArray() )

ace_energy_grid = xss_extractor.extractElectronEnergyGrid()
ace_index = 0
for i in range(0, ace_energy_grid.size ):
    if ace_energy_grid[i] <= energy:
        ace_index = i

energy_0 = ace_energy_grid[ace_index]
energy_1 = ace_energy_grid[ace_index+1]
lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

# -------------------------------------------------------------------------- ##
#  ACE Cutoff Elastic Cross Section
# -------------------------------------------------------------------------- ##
elastic_cross_sections = xss_extractor.extractElasticCrossSection()
cs_0 = elastic_cross_sections[ace_index]
cs_1 = elastic_cross_sections[ace_index+1]

#ace_elastic_cs = cs_0 + (cs_1 - cs_0)*lin_interp
ace_elastic_cs = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)

# -------------------------------------------------------------------------- ##
#  ENDL Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'Au-ENDL' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
endl_data = ENDL.ENDLDataContainer( file_name )

# -------------------------------------------------------------------------- ##
#  ENDL Cutoff Elastic Cross Section
# -------------------------------------------------------------------------- ##
elastic_cs = endl_data.getCutoffElasticCrossSection()
endl_energy_grid = endl_data.getElasticEnergyGrid()

endl_index = 0
for n in range(0, endl_energy_grid.size ):
    if endl_energy_grid[n] <= energy:
        endl_index = n

energy_0 = endl_energy_grid[endl_index]
energy_1 = endl_energy_grid[endl_index+1]
lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

cs_0 = elastic_cs[endl_index]
cs_1 = elastic_cs[endl_index+1]

#endl_elastic_cs = cs_0 + (cs_1 - cs_0)*lin_interp
endl_elastic_cs = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)

# -------------------------------------------------------------------------- ##
#  ENDL Analog Cross Section
# -------------------------------------------------------------------------- ##
analog_cs = endl_data.getTotalElasticCrossSection()

cs_0 = analog_cs[endl_index]
cs_1 = analog_cs[endl_index+1]

#endl_analog_cs = cs_0 + (cs_1 - cs_0)*lin_interp
endl_analog_cs = numpy.exp(numpy.log(cs_0) + numpy.log(cs_1/cs_0)*log_interp)


angular_energy_grid = endl_data.getCutoffElasticAngularEnergyGrid()

angular_index = 0
for n in range(0, angular_energy_grid.size ):
    if angular_energy_grid[n] <= energy:
        angular_index = n

energy_0 = angular_energy_grid[angular_index]
energy_1 = angular_energy_grid[angular_index+1]
lin_interp = ( energy - energy_0 )/( energy_1 - energy_0 )
log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)

pdf_0 = endl_data.getCutoffElasticPDFAtEnergy(energy_0)
cutoff_pdf_0 = pdf_0[0]

pdf_1 = endl_data.getCutoffElasticPDFAtEnergy(energy_1)
cutoff_pdf_1 = pdf_1[0]

cutoff_pdf = cutoff_pdf_0 + (cutoff_pdf_1 - cutoff_pdf_0)*log_interp
cutoff_pdf_lin = cutoff_pdf_0 + (cutoff_pdf_1 - cutoff_pdf_0)*lin_interp
cutoff_pdf_loglog = numpy.exp(numpy.log(cutoff_pdf_0) + numpy.log(cutoff_pdf_1/cutoff_pdf_0)*log_interp)


eta = analog_dist.evaluateMoliereScreeningConstant( energy )

cutoff_cdf = eta/( eta + cutoff_pdf*(eta*1.0e-6 + 1.0e-12) )
cutoff_cdf_lin = eta/( eta + cutoff_pdf*(eta*1.0e-6 + 1.0e-12) )
cutoff_cdf_loglog = eta/( eta + cutoff_pdf*(eta*1.0e-6 + 1.0e-12) )

endl_analog_cs_2 = endl_elastic_cs/cutoff_cdf
endl_analog_cs_3 = endl_elastic_cs/cutoff_cdf_lin
endl_analog_cs_4 = endl_elastic_cs/cutoff_cdf_loglog

# -------------------------------------------------------------------------- ##
#  Elastic Data Comparison
# -------------------------------------------------------------------------- ##

print "\n----- Cutoff Elastic Cross Section at",energy,"-----"

print '\tNative cs = ','%.16e' % native_elastic_cs
print '\tACE cs    = ','%.16e' % ace_elastic_cs
print '\tENDL cs   = ','%.16e' % endl_elastic_cs
rel_diff = (native_elastic_cs-endl_elastic_cs)/endl_elastic_cs
print '\tRel diff (Native vs ENDL) = ','%.16e' % rel_diff
rel_diff = (ace_elastic_cs-endl_elastic_cs)/endl_elastic_cs
print '\tRel diff (ACE vs ENDL)    = ','%.16e' % rel_diff

print "\n----- Analog Elastic Cross Section at",energy,"-----"

print '\tNative cs        = ','%.16e' % native_analog_cs
print '\tENDL cs (tables) = ','%.16e' % endl_analog_cs
print '\tENDL cs (calc)   = ','%.16e' % endl_analog_cs_2
print '\tENDL cs (calc)   = ','%.16e' % endl_analog_cs_3
print '\tENDL cs (calc)   = ','%.16e' % endl_analog_cs_4
rel_diff = (native_analog_cs-endl_analog_cs)/endl_analog_cs
print '\trel diff (Native vs ENDL Tables)    = ','%.16e' % rel_diff
rel_diff = (endl_analog_cs_2-endl_analog_cs)/endl_analog_cs
print '\trel diff (lin-log ENDL Calc vs ENDL Tables) = ','%.16e' % rel_diff
rel_diff = (endl_analog_cs_3-endl_analog_cs)/endl_analog_cs
print '\trel diff (lin-lin ENDL Calc vs ENDL Tables) = ','%.16e' % rel_diff
rel_diff = (endl_analog_cs_4-endl_analog_cs)/endl_analog_cs
print '\trel diff (log-log ENDL Calc vs ENDL Tables) = ','%.16e' % rel_diff



