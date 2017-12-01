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
data_list = cs_list.get( 'H-Native' )

# -------------------------------------------------------------------------- ##
#  Adjoint Elastic Data
# -------------------------------------------------------------------------- ##
adjoint_file_name = datadir + data_list.get( 'adjoint_electroatomic_file_path' )
adjoint_data = Native.AdjointElectronPhotonRelaxationDataContainer( adjoint_file_name )
adjoint_energy_grid = adjoint_data.getAdjointElectronEnergyGrid()


tot_adjoint_elastic_cs = adjoint_data.getAdjointTotalElasticCrossSection()
adjoint_cutoff_cs = adjoint_data.getAdjointCutoffElasticCrossSection()
reduced_cutoff_ratio = adjoint_data.getReducedCutoffCrossSectionRatios()
adjoint_screen_rutherford_cs = adjoint_data.getAdjointScreenedRutherfordElasticCrossSection()
adjoint_screen_rutherford_index = adjoint_data.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex()

adjoint_excitation_cs = adjoint_data.getAdjointAtomicExcitationCrossSection()
adjoint_ionization_cs = adjoint_data.getAdjointElectroionizationCrossSection(1)
adjoint_brem_cs = adjoint_data.getAdjointBremsstrahlungElectronCrossSection()

energies = [ 20.0, 3.2e-2, 1.0e-2]
angles = [-1.0, 0.0, 0.999999, 1.0]

###
###  Coupled Distribution/Reaction Unit Test Check
###

adjoint_dist = Collision.createLinLinLogExactCoupledElasticDistribution(adjoint_data, "Simplified Union", 1e-7)

print "\n------------------------------------------------"
print "\nAdjoint Coupled"
print "------------------------------------------------"

print "\nAdjoint Coupled evaluate"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        pdf = adjoint_dist.evaluate( energy, angle )
        print '\teval(', '%.5e' % angle,') =\t','%.16e' % pdf

print "\nAdjoint Coupled evaluatePDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        pdf = adjoint_dist.evaluatePDF( energy, angle )
        print '\tpdf(', '%.5e' % angle,') =\t','%.16e' % pdf

print "\nAdjoint Coupled evaluateCDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        cdf = adjoint_dist.evaluateCDF( energy, angle )
        print '\tcdf(', '%.5e' % angle,') =\t','%.16e' % cdf


print "\nAdjoint Coupled sample"
print "------------------------------------------------"

random_numbers = [ 0.5, 0.5 ]
for energy in energies:
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
    print "\n---- energy =", energy, "----"
    e_out, mu = adjoint_dist.sample( energy )
    print '\tmu =\t','%.16e' % mu
    print '\te_out =\t','%.16e' % e_out

###
###  Cutoff Distribution/Reaction Unit Test Check
###
print "\n------------------------------------------------"
print "Adjoint Cutoff"
print "------------------------------------------------"

angles = [-1.0, 0.0, 0.9, 0.9001]
adjoint_dist = Collision.createLinLinLogExactCutoffElasticDistribution(adjoint_data, 0.9, 1e-7)

print "\nAdjoint Cutoff evaluate"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        pdf = adjoint_dist.evaluate( energy, angle )
        print '\teval(', '%.5e' % angle,') =\t','%.16e' % pdf

print "\nAdjoint Cutoff evaluatePDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        pdf = adjoint_dist.evaluatePDF( energy, angle )
        print '\tpdf(', '%.5e' % angle,') =\t','%.16e' % pdf

print "\nAdjoint Cutoff evaluateCDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        cdf = adjoint_dist.evaluateCDF( energy, angle )
        print '\tcdf(', '%.5e' % angle,') =\t','%.16e' % cdf


print "\nAdjoint Cutoff sample"
print "------------------------------------------------"

random_numbers = [ 0.5, 0.5 ]
for energy in energies:
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
    print "\n---- energy =", energy, "----"
    e_out, mu = adjoint_dist.sample( energy )
    print '\tmu =\t','%.16e' % mu
    print '\te_out =\t','%.16e' % e_out

###
###  Hybrid Distribution/Reaction Unit Test Check
###
print "\n------------------------------------------------"
print "Adjoint Hybrid"
print "------------------------------------------------"

energies = [ 3.2e-2, 1.0e-2]
angles = [0.0, 0.9, 0.9001]
adjoint_dist = Collision.createLinLinLogExactHybridElasticDistribution(adjoint_data, 0.9, 1e-7)

print "\nAdjoint Hybrid evaluate"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        pdf = adjoint_dist.evaluate( energy, angle )
        print '\teval(', '%.5e' % angle,') =\t','%.16e' % pdf

print "\nAdjoint Hybrid evaluatePDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        pdf = adjoint_dist.evaluatePDF( energy, angle )
        print '\tpdf(', '%.5e' % angle,') =\t','%.16e' % pdf

print "\nAdjoint Hybrid evaluateCDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        cdf = adjoint_dist.evaluateCDF( energy, angle )
        print '\tcdf(', '%.5e' % angle,') =\t','%.16e' % cdf


print "\nAdjoint Hybrid sample"
print "------------------------------------------------"

random_numbers = [ 0.5, 0.5 ]
for energy in energies:
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
    print "\n---- energy =", energy, "----"
    e_out, mu = adjoint_dist.sample( energy )
    print '\tmu =\t','%.16e' % mu


###
###  Moment Preserving Distribution/Reaction Unit Test Check
###
print "\n------------------------------------------------"
print "Adjoint Moment Preserving"
print "------------------------------------------------"

angular_energy_grid = adjoint_data.getAdjointElasticAngularEnergyGrid()

energies = [ 3.2e-2, 1.0e-2]
angles = [0.0, 0.9, 0.9001]
adjoint_dist = Collision.createLinLinLogExactMomentPreservingElasticDistribution(adjoint_data, 0.9, 1e-7)

discrete_angles = adjoint_data.getAdjointMomentPreservingElasticDiscreteAngles( 3.2e-2 )
lower_discrete_angles = adjoint_data.getAdjointMomentPreservingElasticDiscreteAngles( 8e-3 )
upper_discrete_angles = adjoint_data.getAdjointMomentPreservingElasticDiscreteAngles( 1.6e-2 )

print '%.16e' % lower_discrete_angles[0]
print '%.16e' % lower_discrete_angles[1]
print '%.16e' % upper_discrete_angles[0]
print '%.16e' % upper_discrete_angles[1]
print '%.16e' % discrete_angles[0]
print '%.16e' % discrete_angles[1]

angles = [ lower_discrete_angles[0], lower_discrete_angles[1],upper_discrete_angles[0],upper_discrete_angles[1],discrete_angles[0],discrete_angles[1] ]
print "\nAdjoint Moment Preserving evaluate"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        pdf = adjoint_dist.evaluate( energy, angle )
        print '\teval(', '%.16e' % angle,') =\t','%.16e' % pdf

print "\nAdjoint Moment Preserving evaluatePDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        pdf = adjoint_dist.evaluatePDF( energy, angle )
        print '\tpdf(', '%.16e' % angle,') =\t','%.16e' % pdf

print "\nAdjoint Moment Preserving evaluateCDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for angle in angles:
        cdf = adjoint_dist.evaluateCDF( energy, angle )
        print '\tcdf(', '%.16e' % angle,') =\t','%.16e' % cdf


print "\nAdjoint Moment Preserving sample"
print "------------------------------------------------"

random_numbers = [ 0.0, 1.0-1e-15, 0.0, 1.0-1e-15 ]
for energy in energies:
    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
    print "\n---- energy =", energy, "----"
    e_out, mu = adjoint_dist.sample( energy )
    print '\tmu =\t','%.16e' % mu
    e_out, mu = adjoint_dist.sample( energy )
    print '\tmu =\t','%.16e' % mu

###
###  Bremsstrahlung Distribution/Reaction Unit Test Check
###
print "\n------------------------------------------------"
print "Adjoint Bremsstrahlung"
print "------------------------------------------------"

energies = [ 1e-5, 1.0e-3]
e_outs = [0.0, 1.0, 20.0]
adjoint_dist = Collision.createLinLinLogCorrelatedBremsstrahlungDistribution(adjoint_data, 1e-7)

print "\nAdjoint Bremsstrahlung evaluate"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for e_out in e_outs:
        if e_out < energy:
            e_out = 2.0*energy
        pdf = adjoint_dist.evaluate( energy, e_out )
        print '\teval(', '%.16e' % e_out,') =\t','%.16e' % pdf

print "\nAdjoint Bremsstrahlung evaluatePDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for e_out in e_outs:
        if e_out < energy:
            e_out = 2.0*energy
        pdf = adjoint_dist.evaluatePDF( energy, e_out )
        print '\teval(', '%.16e' % e_out,') =\t','%.16e' % pdf

print "\nAdjoint Bremsstrahlung evaluateCDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for e_out in e_outs:
        if e_out < energy:
            e_out = 2.0*energy
        pdf = adjoint_dist.evaluateCDF( energy, e_out )
        print '\teval(', '%.16e' % e_out,') =\t','%.16e' % pdf


print "\nAdjoint Bremsstrahlung sample"
print "------------------------------------------------"

energy = 1e-5
print "\n---- energy =", energy, "----"

random_numbers = [ 0.0, 0.0, 0.0 ]
Prng.RandomNumberGenerator.setFakeStream(random_numbers)
e_out, mu = adjoint_dist.sample( energy )
print '\te_out =\t','%.16e' % e_out

random_numbers = [ 1.0-1e-15, 1.0-1e-15, 0.0 ]
Prng.RandomNumberGenerator.setFakeStream(random_numbers)
e_out, mu = adjoint_dist.sample( energy )
print '\te_out =\t','%.16e' % e_out

energy = 1e-4
print "\n---- energy =", energy, "----"

random_numbers = [ 0.5, 0.5, 1.0-1e-15, 0.0 ]
Prng.RandomNumberGenerator.setFakeStream(random_numbers)
e_out, mu = adjoint_dist.sample( energy )
print '\te_out =\t','%.16e' % e_out


###
###  Electroionization Distribution/Reaction Unit Test Check
###
print "\n------------------------------------------------"
print "Adjoint Electroionization"
print "------------------------------------------------"

energies = [ 1e-5, 1.0e-3]
e_outs = [0.0, 1.0, 20.0]

shell = 1
binding_energy = adjoint_data.getSubshellBindingEnergy( shell )
adjoint_dist = Collision.createLinLinLogCorrelatedElectroionizationSubshellDistribution(adjoint_data, shell, binding_energy, 1e-7)

print 'Binding energy =\t','%.16e' % binding_energy

print "\nAdjoint Electroionization evaluate"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for e_out in e_outs:
        if e_out < energy:
            e_out = 2.0*energy + binding_energy
        pdf = adjoint_dist.evaluate( energy, e_out )
        print '\teval(', '%.16e' % e_out,') =\t','%.16e' % pdf

print "\nAdjoint Electroionization evaluatePDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for e_out in e_outs:
        if e_out < energy:
            e_out = 2.0*energy + binding_energy
        pdf = adjoint_dist.evaluatePDF( energy, e_out )
        print '\teval(', '%.16e' % e_out,') =\t','%.16e' % pdf

print "\nAdjoint Electroionization evaluateCDF"
print "------------------------------------------------"

for energy in energies:
    print "\n---- energy =", energy, "----"
    for e_out in e_outs:
        if e_out < energy:
            e_out = 2.0*energy + binding_energy
        pdf = adjoint_dist.evaluateCDF( energy, e_out )
        print '\teval(', '%.16e' % e_out,') =\t','%.16e' % pdf


print "\nAdjoint Electroionization sample"
print "------------------------------------------------"

random_numbers = [ 0.0, 1.0-1e-15 ]
for energy in energies:
    print "\n---- energy =", energy, "----"

    Prng.RandomNumberGenerator.setFakeStream(random_numbers)
    e_out, mu = adjoint_dist.sample( energy )
    print '\te_out =\t','%.16e' % e_out, '\tmu =\t','%.16e' % mu
    e_out, mu = adjoint_dist.sample( energy )
    print '\te_out =\t','%.16e' % e_out, '\tmu =\t','%.16e' % mu


