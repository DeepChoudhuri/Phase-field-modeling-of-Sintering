import numpy as np
import matplotlib.pyplot as plt

# For plotting curves using NINE particles
#f = np.loadtxt('time-step-energy.txt')

#f = np.loadtxt('type-6_Case-I/A-16_B-1_kp-0.3_kn-6.7/small_circle_data_structureType-6.txt')
#f = np.loadtxt('type-6_Case-I/A-16_B-1_kp-13.7_kn-2/small_circle_data_structureType-6.txt')
#f = np.loadtxt('type-6_Case-I/A-16_B-1_kp-2.5_kn-4.5/small_circle_data_structureType-6.txt')

#f = np.loadtxt('type-8_case-IVb/A-16_B-1_kp-0.3_kn-6.7/small_circle_data_structureType-8.txt')
#f = np.loadtxt('type-8_case-IVb/A-16_B-1_kp-13.7_kn-2/small_circle_data_structureType-8.txt')
#f = np.loadtxt('type-8_case-IVb/A-16_B-1_kp-2.5_kn-4.5/small_circle_data_structureType-8.txt')

#f = np.loadtxt('type-7_case-III/A-16_B-1_kp-2.5_kn-4.5/small_circle_data_structureType-7.txt')

#f = np.loadtxt('type-2_Case-II/A-16_B-1_kp-0.3_kn-6.7/small_circle_data_structureType-2.txt')
#f = np.loadtxt('type-2_Case-II/A-16_B-1_kp-13.7_kn-2/small_circle_data_structureType-2.txt')

#f = np.loadtxt('type-5_Case-V/A-16_B-1_kp-13.7_kn-2/small_circle_data_structureType-5.txt')

#f = np.loadtxt('type-9_Case-VI/A-16_B-1_kp-13.7_kn-2/small_circle_data_structureType-9.txt')

f = np.loadtxt('small_circle_data_structureType-8.txt')

#plt.ylable("Free Energy")
delX = 0.4
plt.figure(1)
plt.subplot(211).set_ylabel('Total neck length')
#plt.plot(f[:,1],np.sqrt(2*delX**2)*f[:,4],'b-')
plt.plot(f[:,1],np.sqrt(2*delX**2)*f[:,4],'bo')
print("Maximum neck length ", np.max(np.sqrt(2*delX**2)*f[:,4]))

plt.subplot(212).set_ylabel('Particle area')
plt.plot(f[:,1],f[:,3],'ro')

plt.xlabel('No. iterations')

plt.show()
