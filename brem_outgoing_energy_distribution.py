#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import numpy as numpy

Utility.initFrensiePrng()

#datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/research/frensie-repos/lkersting/src/packages/test_files/'

# -------------------------------------------------------------------------- ##
#  Electron Data
# -------------------------------------------------------------------------- ##

file_name = datadir + 'native/test_epr_1_native.xml'
native_data = Native.ElectronPhotonRelaxationDataContainer(file_name)

energy_grid = list (native_data.getBremsstrahlungEnergyGrid() )

for i in range(len(energy_grid)):
  print "\n", i, "\n"
  energy = energy_grid[i]
  photon_energy = list( native_data.getBremsstrahlungPhotonEnergy( energy ) )
  photon_pdf = list( native_data.getBremsstrahlungPhotonPDF( energy ) )

  primary_energy = photon_energy[::-1]
  primary_pdf = photon_pdf[::-1]


  # Evaluate the energy of the primary electron energy at the incoming energy
  primary_energy = map(lambda x: energy - x, primary_energy )
  primary_energy[0] += 1e-10

  print energy
  print photon_energy[-1], photon_energy[-2], photon_energy[0]
  print energy - photon_energy[-1], energy - photon_energy[-2], energy - photon_energy[0]
  print primary_energy[0], primary_energy[1], primary_energy[-1]