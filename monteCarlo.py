#Import modules
import numpy as np
import matplotlib.pyplot as plt
import itertools


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

    # Get system potential energy
    def getHam(self):
        ham = 0
        for i in range(0, len(self.particles) - 1):
            for j in range(i + 1, len(self.particles)):
                r = np.linalg.norm(self.particles[j].pos - self.particles[i].pos)
                if r < cutoff:
                    ham += lenJones(r)
        return ham
    
    # Display particles in system
    def displayParticles(self, particles=None):
        if particles is None:
            particles = self.particles
    
        xList = [p.pos[0] for p in particles]
        yList = [p.pos[1] for p in particles]
        zList = [p.pos[2] for p in particles]
        
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(xList, yList, zList, marker="o")
        plt.show()
    
    # Get particles within a specified region
    
    def get_particles_in_region(self, coords1, coords2):
        selected_particles = []
    
        for particle in self.particles:
            if all(min_coords <= particle.pos) and all(particle.pos <= max_coords):
                selected_particles.append(particle)
        
        return selected_particles

class particle():
	
	def __init__(self):
		self.pos = np.array([np.random.uniform(0,cubeSide),np.random.uniform(0,cubeSide),np.random.uniform(0,cubeSide)])

B = box()
print(lenJones(10))
print(B.getHam())


# generate boxes 27 (3-dimensions) 

# determining potentals, pair potentials 

# Define region boundries

R_c = 50
min_coords = np.array([0, 0, 0]) 
inital_array = [R_c, cubeSide, cubeSide]
volumes = []
region_box = B

for _ in range(0, 3):
    max_coords = np.array(np.roll((inital_array), _))
    particles_in_region = region_box.get_particles_in_region(min_coords, max_coords)
    volumes.append(particles_in_region)

max_coords = np.array([cubeSide, cubeSide, cubeSide])
inital_array = [(cubeSide-R_c), 0, 0]

for _ in range(0,3):
    min_coords = np.array(np.roll((inital_array), _))
    particles_in_region = region_box.get_particles_in_region(min_coords, max_coords)
    volumes.append(particles_in_region)

for i in range(len(volumes)):
    print("Number of particles in region:" + str(i) + " " + str(len(volumes[i])))

B.displayParticles(volumes[0])

# calculate boundry potentals, bend potential
# obtain hamitoian new - old
# determine probabilty 
# compare probabilty with random number