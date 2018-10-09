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

# Get the atomic shells
subshells = native_data.getSubshells()
shell = subshells[0]
binding_energy = native_data.getSubshellBindingEnergy( shell )

energy_grid = list (native_data.getElectroionizationEnergyGrid(shell) )

for i in range(len(energy_grid)):
  print "\n", i, "\n"
  energy = energy_grid[i]
  knock_on_energy = list( native_data.getElectroionizationRecoilEnergy( shell, energy ) )
  knock_on_pdf = list( native_data.getElectroionizationRecoilPDF( shell, energy ) )

  primary_energy = knock_on_energy[:-1][::-1]
  primary_pdf = knock_on_pdf[:-1][::-1]
  print primary_energy


  # If the grid point is not greater than the binding energy then fix it
  if energy_grid[i] <= binding_energy:
    # Use the max knock-on energy to calculate the incoming energy
    energy_2k = 2.0*knock_on_energy[-1]
  if energy > binding_energy:
    energy_2k = energy - binding_energy
    knock_on_energy[-1] = energy_2k/2.0

  # Evaluate the energy of the primary electron energy at the incoming energy
  primary_energy = map(lambda x: energy_2k - x, primary_energy )
  print primary_energy

  print "diff = ", knock_on_energy[-1]*2, knock_on_energy[-2] + primary_energy[0], knock_on_energy[-1]*2 -( knock_on_energy[-2] + primary_energy[0])

  energy_bins = knock_on_energy + primary_energy
  pdfs = knock_on_pdf + primary_pdf

  print energy_bins
  size = len(primary_energy)
  for i in range(size):
    print "energy add =", energy_bins[i] + energy_bins[-i-1], "\tpdf diff =\t", pdfs[i] - pdfs[-i-1]



  # # If the grid point is not greater than the binding energy then fix it
  # if energy_grid[i] <= binding_energy:
  #   # Use the max knock-on energy to calculate the incoming energy
  #   energy = 2.0*knock_on_energy[-1] + binding_energy
  #   energy_grid[i] = energy
  # if energy > binding_energy:
  #   print binding_energy - energy
  #   knock_on_energy[-1] = (energy - binding_energy)/2.0

  # # Evaluate the energy of the primary electron energy at the incoming energy
  # primary_energy = map(lambda x: energy - x - binding_energy, primary_energy )
  # print primary_energy

  # print "diff = ", knock_on_energy[-1]*2, knock_on_energy[-2] + primary_energy[0], knock_on_energy[-1]*2 -( knock_on_energy[-2] + primary_energy[0])

  # energy_bins = knock_on_energy + primary_energy
  # pdfs = knock_on_pdf + primary_pdf

  # print energy_bins
  # size = len(primary_energy)
  # for i in range(size):
  #   print "energy add =", energy_bins[i] + energy_bins[-i-1], "\tpdf diff =\t", pdfs[i] - pdfs[-i-1]