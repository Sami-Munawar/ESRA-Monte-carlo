#Import modules
import numpy as np
import matplotlib.pyplot as plt

#Units: angstroms, daltons, kelvin, mole
#Boltzmann constant defined as 1

#Define physics constants, number of particles in simulation, number of moles and simulation parameters
avogadro = 6.02214076e+23
R = 8.31446261815324
N = 500
moles = N/avogadro
T = 298
#Molar volume in SI units
molVol = 0.024789570296023
#Convert volume into a cube side in angstroms cubed
cubeSide = (moles*molVol/1e-30)**(1/3)

#Units in simulation: angstroms and daltons expressed in SI units
angstrom = 1e-10
dalton = 1.66053906892e-27

#Define Lennard-Jones potential constants for Xenon
#Boltzmann constant, in J
kB = 1.380649e-23
epsilon = 229
sigma = 4.06
sigma_12 = sigma**12
sigma_6 = sigma**6
#Cutoff radius
cutoff = 5*sigma
alpha = 4*epsilon*(12*(sigma_12)/(cutoff**13)-6*(sigma_6)/(cutoff**7))
beta = -4*epsilon*((sigma_12)/(cutoff**12)-(sigma_6)/(cutoff**6))-alpha*cutoff
#The maximum amount a particle can be displaced in 1 step
maxDis = 3*sigma

#Helper methods

#Get probability of system changing state
def getProb(initHam,newHam):
	return np.exp(-(newHam-initHam)/(T))

#Lennard-Jones potential (truncated)
def lenJones(radius):
	energy = 4*epsilon*((sigma_12)/(radius**12)-(sigma_6)/(radius**6))+alpha*radius+beta
	return energy

#Generate random displacement for a particle
def getDis():
	disMag = np.random.uniform(0,maxDis)
	while True:
		x = np.random.normal(0,1)
		y = np.random.normal(0,1)
		z = np.random.normal(0,1)
		if (x != 0 and y != 0) and z != 0:
			break
	normFactor = (x**2+y**2+z**2)**(1/2)
	return (disMag/normFactor)*np.array([x,y,z])

#Classes
class box():
	
	def __init__(self):
		self.particles = np.array([particle() for i in range(N)])
	
	#Get system potential energy
	def getHam(self):
		ham = 0
		for i in range(0,len(self.particles)-1):
			for j in range(i+1,len(self.particles)):
				r = np.linalg.norm(self.particles[j].pos-self.particles[i].pos)
				if r < cutoff:
					ham += lenJones(r)
		return ham
	
	#Display particles in system
	def displayParticles(self):
		xList = [particle.pos[0] for particle in self.particles]
		yList = [particle.pos[1] for particle in self.particles]
		zList = [particle.pos[2] for particle in self.particles]
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')
		ax.scatter(xList,yList,zList,marker="o")
		plt.show()

class particle():
	
	def __init__(self):
		self.pos = np.array([np.random.uniform(0,cubeSide),np.random.uniform(0,cubeSide),np.random.uniform(0,cubeSide)])

B = box()
print(lenJones(10))
print(B.getHam())
B.displayParticles()

# generate boxes 27 (3-dimensions) 

boxes = [] 
for box in range(27): # adds the 27 box objects to boxes list
	boxes.append(box())

# define pos of box objects




# copy the box objects



# 
# determining potentals, pair potentials 
# calculate boundry potentals, bend potential
# obtain hamitoian new - old
# determine probabilty 
# compare probabilty with random number