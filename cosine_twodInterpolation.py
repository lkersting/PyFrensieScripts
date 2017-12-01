#! /usr/bin/env python
import numpy as np

print "\n\t### Test Answers ###"

# Values
x0 = 0.1; x1 = 1.0; y = 0.9

y_min_0 = 0.5; y_min_1 = 0.999999

z_min_0 = 1.0; z_min_1 = 10.0

y_max_0 = 0.0; y_max_1 = 0.999999

z_max_0 = 5.0; z_max_1 = 0.5

# Processed values
px_0 = np.log(x0); px_1 = np.log(x1); py = np.log(1.0 - y)
# print "px_0 =",px_0,"\tpx_1 =",px_1
# print "py =",py
py_min_0 = np.log(1.0 - y_min_0); py_min_1 = np.log(1.0 - y_min_1)
#print "py_min_0 =",py_min_0,"\tpy_min_1 =",py_min_1
pz_min_0 = np.log(z_min_0); pz_min_1 = np.log(z_min_1)
# print "pz_min_0 =",pz_min_0,"\tpz_min_1 =",pz_min_1
py_max_0 = np.log(1.0 - y_max_0); py_max_1 = np.log(1.0 - y_max_1)

pz_max_0 = np.log(z_max_0); pz_max_1 = np.log(z_max_1)

# Interpolate the processed min and max z
pz_0 = pz_min_0 + ( pz_min_1 - pz_min_0 )*( py - py_min_0 )/( py_min_1 - py_min_0 )
pz_1 = pz_max_0 + ( pz_max_1 - pz_max_0 )*( py - py_max_0 )/( py_max_1 - py_max_0 )

# Unprocessed min and max z
z_0 = np.exp( pz_0 ); z_1 = np.exp( pz_1 )
  
# Interpolate the z value
processed_slope = ( pz_1 - pz_0 )/(px_1 - px_0)
# print "processed_slope = ",processed_slope
print "\n--- interpolate_separate_tuple_grids ---"

x_values = [0.3, 0.1, 1.0]

for x in x_values:
    #process x value
    px = np.log(x)
    print "px = ",px
    # Interpolate the z value
    pz = pz_0 + (px - px_0)*processed_slope

    # Unprocessed z
    z = np.exp( pz )
    print "\tx =", x, "\tz =",'%.16e'% z


print "\n--- interpolateUnitBase_separate_tuple_grids ---"
true_end = 6.9314718055994529e-01
max_y = 6.9314718055994540e-01
tol = 1e-3
fuzzy = true_end*(1+tol)
print "max_y = ",max_y
print "fuzzy = ",fuzzy
py_max = np.log( 1.0 - (1.0-1e-15) )
py_min = np.log( 2.0 )
print py_min - py_max
# x values
x0 = 0.1; x1 = 1.0; x_values = [0.3, 0.1, 1.0]
x_values = [0.3]
# y values
y_values = [ 0.9, -1.0, 0.999999 ]
y_values = [0.9, -1.0]
y_min_0_values = [ 0.5, -1.0, 0.5 ]
y_min_1_values = [ 0.999999, 0.0, 0.999999 ]
y_max_0_values = [ 0.0, -1.0, 0.0 ]
y_max_1_values = [ 0.999999, 0.0, 0.999999 ]

# z values
z_min_0_values = [ 0.1, 1e-3, 0.1 ]
z_min_1_values = [ 1.0, 1e-2, 1.0 ]
z_max_0_values = [ 0.1, 1e-2, 0.1 ]
z_max_1_values = [ 1.0, 1e-1, 1.0 ]

for j in range(len(x_values)):
  # Set x value
  x = x_values[j]

  # Processed x value
  px = np.log(x)

  print "\nx =", x
  for i in range(len(y_values)):
    # Set vlaues
    y = y_values[i]
    y_min_0 = y_min_0_values[i]; y_min_1 = y_min_1_values[i]
    y_max_0 = y_max_0_values[i]; y_max_1 = y_max_1_values[i]

    z_min_0 = z_min_0_values[i]; z_min_1 = z_min_1_values[i]
    z_max_0 = z_max_0_values[i]; z_max_1 = z_max_1_values[i]

    # Processed values
    print "y =",y
    py = np.log(1.0 - y)
    print "py =", py
    py_min_0 = np.log(1.0 - y_min_0); py_min_1 = np.log(1.0 - y_min_1)
    py_max_0 = np.log(1.0 - y_max_0); py_max_1 = np.log(1.0 - y_max_1)

    pz_min_0 = np.log(z_min_0); pz_min_1 = np.log(z_min_1)
    pz_max_0 = np.log(z_max_0); pz_max_1 = np.log(z_max_1)

    # Get grid lengths
    y_min = -1.0; y_max = 0.999999
    py_min = np.log( 1.0 - y_min ); py_max = np.log( 1.0 - y_max )
    print "py_min =", py_min
    L = -(py_max - py_min)
    print "L =",L
    # Get eta
    eta = (py - py_min)/L
    print "eta =", eta
    print 1.0+eta
    print (py - py_max)/L
    # Get y values
    py0_min = py_min + L*eta; py1_min = py0_min
    print "py0_min =", py0_min
    print py_max + L*(py - py_max)/L
    # Interpolate the processed min and max z
    pz_0 = pz_min_0 + ( pz_min_1 - pz_min_0 )*( py - py_min_0 )/( py_min_1 - py_min_0 )
    pz_1 = pz_max_0 + ( pz_max_1 - pz_max_0 )*( py - py_max_0 )/( py_max_1 - py_max_0 )

    # Scale min adn max
    scaled_pz_0 = pz_0*L; scaled_pz_1 = pz_1*L

    # Interpolate the z value
    processed_slope = ( scaled_pz_1 - scaled_pz_0 )/(px_1 - px_0)
    print "processed_slope =", processed_slope
    # Interpolate the z value
    pz = ( scaled_pz_0 + (px - px_0)*processed_slope )/L

    # Unprocessed z
    z = np.exp( pz )
    print "\ty =", y, "\tz =",'%.16e'% z