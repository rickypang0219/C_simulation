import numpy as np 
import scipy.linalg as la 
import matplotlib.pyplot as plt 
import pandas as pd
import time as time

'''
Objective: 
0. Investigate the Markov Chain Monte Carlo
1. Compute the magnetization curve VS inverse tempreture  
2. Compute the auto-correlation and see whether a certain Markov chain can produce less correlated sample 
3. After finishing this, try to convert in C
'''



# Create a class of Ising Model
# In this ising class, it should include
# 1. Parameter --> Size , Boundary condidtion, strength parameter J, \mu 
# 2. From Strengths --> Construct 2D Ising Hamiltonian ( No need to construct, only need to find out the energy difference of spin flip of a site)
# 3. MCMC generates different samples of configguration in Ising model
# 3.1 what are the stationary distribution 
# 3.2 Suppose at temperature T, I sample 1000 times configuration, is it enough ? 
# 3.3 What are we sampling? Sample the spin configurations? 

class Ising: 

    def __init__(self, size, J, mu, Temp_max=6): 
        self.size = int(size) 
        self.J = J 
        self.mu = mu 
        self.Temp_max = Temp_max
        self.lattice = np.random.choice([-1, 1], size=(self.size, self.size), p=[0.5, 0.5])
        self.dT = 0.2
        self.sweep = 1000
        self.bins = 1
        self.T_array = np.arange(1E-1, self.Temp_max, self.dT)
        self.Magnetization_Matrix = np.zeros((len(self.T_array), self.bins * self.sweep))
        self.Nbins = 5
        self.Nsweep = 1000
        self.num_T = len(self.T_array)
        self.all_M=np.zeros([self.num_T,self.Nbins*self.Nsweep])
        self.aver_M=np.zeros(self.num_T)
    

    def get_energy_diff(self, i:int, j:int) -> float:
        current_site = self.lattice[i,j]
        adjacent_site = self.lattice[(i+1)%self.size, j] + self.lattice[(i-1)%self.size,j] + self.lattice[i, (j+1)%self.size] + self.lattice[i, (j-1)%self.size]
        energy_diff = 2 * current_site * ( self.J * adjacent_site + self.mu )
        return energy_diff

   
    def get_Metropolis_sampling(self):
        '''
        Using the energy difference, compute the acceptance ratio 
        Return configuration
        '''
        dt = 0
        for T in self.T_array:  
            for sweep in range(self.bins*self.sweep):
                for _ in range(self.size**2):
                    i = np.random.randint(0,self.size)
                    j = np.random.randint(0,self.size)
                    delta_E = self.get_energy_diff(i,j)
                    if np.exp(- delta_E/T) > np.random.random():
                        self.lattice[i,j] = -self.lattice[i,j]
    
                self.Magnetization_Matrix[dt, sweep] = np.sum(self.lattice)/(self.size**2)
            dt += 1
            print("Complete a sweep at time:", T)

    # def get_Metropolis_sampling(self): 
    #     nT = 0
    #     for T in self.T_array:
    #         nt2=0
    #         start_time = time.time()
    #         for nb in range(self.Nbins):        
    #             for ns in range(self.Nsweep):
    #                 for nl in range(self.size*self.size):
    #                     nx=np.random.randint(self.size)
    #                     ny=np.random.randint(self.size)
    #                     delta_E= self.get_energy_diff(nx,ny)
    #                     if np.exp(-delta_E/T)>np.random.random():
    #                         self.lattice[nx,ny] = - self.lattice[nx,ny]
    #                        
    #                 self.all_M[nT,nt2]=sum(sum(self.lattice))/self.size/self.size
    #                 nt2 +=1
    #                
    #         self.aver_M[nT]=sum(self.all_M[nT,:])/self.Nbins/self.Nsweep
    #         nT+=1
    #         print("--- %s seconds ---" % (time.time() - start_time))
    #
    
    def get_magnetization(self):
        self.get_Metropolis_sampling()
        return np.mean(self.Magnetization_Matrix, axis=1)
        # return self.aver_M






if __name__ == "__main__":
    model = Ising(16,1,0,5)
    # print(model.lattice)
    # print(model.calculate_energy())
    m_avg  = model.get_magnetization()
    # print(model.lattice)
    # df = pd.DataFrame({"Magnetization":m_avg, "Temp":model.T_array})
    # df.to_csv("output.csv")
    plt.plot(model.T_array, np.abs(m_avg), "*")
    plt.show()
