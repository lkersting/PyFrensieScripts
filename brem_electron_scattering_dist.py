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
#  Brem Data
# -------------------------------------------------------------------------- ##
data_list = cs_list.get( 'Pb-Native' )
file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )


energy_grid = native_data.getBremsstrahlungEnergyGrid()
linlinlin_brem_dist = Collision.createBremsstrahlungDistribution(native_data, "LinLinLin", True, False, 1e-7)
linlinlog_brem_dist = Collision.createBremsstrahlungDistribution(native_data, "LinLinLog", True, True, 1e-7)
logloglog_brem_dist = Collision.createBremsstrahlungDistribution(native_data, "LogLogLog", True, False, 1e-7)

energy = 0.0009

index = 0
for i in range(0, energy_grid.size ):
    if energy_grid[i] <= energy:
        index = i

energy_0 = energy_grid[index]
energy_1 = energy_grid[index+1]
print energy_0
print native_data.getBremsstrahlungPhotonEnergy( energy_0 )
print energy_1
print native_data.getBremsstrahlungPhotonEnergy( energy_1 )

cdf_value = [0.5]
for cdf in cdf_value:
    random_numbers = [ cdf, cdf, cdf, cdf, cdf, cdf, cdf ]
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)

    print "\n\t--- Lower Incoming Electron Energy ",energy_0," ---"
    outgoing_E_0, angle_0 = linlinlin_brem_dist.sample( energy_0 )
    linlinlog_outgoing_E_0, linlinlog_angle_0 = linlinlog_brem_dist.sample( energy_0 )
    print "\tLin-Lin-Lin: E_0 = ",'%.18e' % outgoing_E_0,"\tangle_0 = ",'%.18e' % angle_0
    print "\tLin-Lin-Log: E_0 = ",'%.18e' % linlinlog_outgoing_E_0,"\tangle_0 = ",'%.18e' % linlinlog_angle_0

    print "\n\t--- Upper Incoming Electron Energy ",energy_1," ---"
    outgoing_E_1, angle_1 = linlinlin_brem_dist.sample( energy_1 )
    linlinlog_outgoing_E_1, linlinlog_angle_1 = linlinlog_brem_dist.sample( energy_1 )
    print "\tLin-Lin-Lin: E_1 = ",'%.18e' % outgoing_E_1,"\tangle_1 = ",'%.18e' % angle_1
    print "\tLin-Lin-Log: E_1 = ",'%.18e' % linlinlog_outgoing_E_1,"\tangle_1 = ",'%.18e' % linlinlog_angle_1


    E_in_log_interp = numpy.log(energy/energy_0)/numpy.log(energy_1/energy_0)
    E_in_lin_interp = (energy-energy_0)/(energy_1-energy_0)

    loglog_outgoing_E = numpy.exp(numpy.log(outgoing_E_0) + numpy.log(outgoing_E_1/outgoing_E_0)*E_in_log_interp)
    linlog_outgoing_E = outgoing_E_0 + (outgoing_E_1 - outgoing_E_0)*E_in_log_interp
    linlin_outgoing_E = outgoing_E_0 + (outgoing_E_1 - outgoing_E_0)*E_in_lin_interp

    print "\n\t--- Incoming Electron Energy ",energy," ---"
    print "\tLinLin Interp Outgoing Energy: ",'%.18e' % linlin_outgoing_E
    print "\tLinLog Interp Outgoing Energy: ",'%.18e' % linlog_outgoing_E
    print "\tLogLog Interp Outgoing Energy: ",'%.18e' % loglog_outgoing_E

    sample_linlinlin_E, sample_linlinlin_angle = linlinlin_brem_dist.sample( energy )
    print "\tLin-Lin-Lin: Outgoing Energy: = ",'%.18e' % sample_linlinlin_E,"\tangle = ",'%.18e' % sample_linlinlin_angle

    sample_linlinlog_E, sample_linlinlog_angle = linlinlog_brem_dist.sample( energy )
    print "\tLin-Lin-Log: Outgoing Energy: = ",'%.18e' % sample_linlinlog_E,"\tangle = ",'%.18e' % sample_linlinlog_angle

    sample_logloglog_E, sample_logloglog_angle = logloglog_brem_dist.sample( energy )
    print "\tLog-Log-Log: Outgoing Energy: = ",'%.18e' % sample_logloglog_E,"\tangle = ",'%.18e' % sample_logloglog_angle



