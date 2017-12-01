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
data_list = cs_list.get( 'Al-Native' )

### -------------------------------------------------------------------------- ##
###  Forward Elastic Unit Test Check
### -------------------------------------------------------------------------- ##
native_file_name = datadir + data_list.get( 'electroatomic_file_path' )
native_data = Native.ElectronPhotonRelaxationDataContainer( native_file_name )


###
###  Analog vs Cutoff DCS
###
cutoff_reaction = Collision.createCutoffElasticReaction(native_data, 1.0, True, True, 1e-15)
analog_reaction = Collision.createAnalogElasticReaction(native_data, True, True, 1e-15)


energy = 15.7
cosines = numpy.linspace(-1.0,0.999999,100)
for cosine in cosines:
    analog_dcs = analog_reaction.getDifferentialCrossSection(1e-5, cosine)
    cutoff_dcs = cutoff_reaction.getDifferentialCrossSection(1e-5, cosine)
    diff = analog_dcs - cutoff_dcs
    print "difference (",cosine,") =", diff

