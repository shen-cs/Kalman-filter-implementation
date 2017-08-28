import numpy as np
from KFilter import KFilter

measurements = np.ones(shape=(5, 1)) # measurements of accel
initial_x = [0]
dt = 1
z = np.array([[0]])
x = np.array([[0.],  [0.], [0.]]) # initial state (acceleration, velocity, position)
u = np.array([[0.], [0.], [0.]]) 
P =  np.array([[0, 0, 0], [0, 1000, 0], [0, 0, 1000]])
F =  np.array([[1., dt, 0.5*dt**2.], [0, 1, dt], [0, 0, 1]]) 
H =  np.array([[0., 0., 1.]])
R =  np.array([[10.]]) 
I =  np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) # 4d identity matrix

kfilter = KFilter(x, P, F, u, z, H, R)

for i in range(len(measurements)):
  kfilter.filterOnce(measurements[i].reshape(1, 1))
print 'x= '
print kfilter.x
print 'P= '
print kfilter.P