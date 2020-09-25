#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that simulates and outputs individual gravitational waveforms. Partially
based on the work of Huerta et al. (2017) and Buskirk et al. (2019), as well as
scripts by JJ Eldridge.
"""

import numpy as np
import matplotlib.pyplot as plt

"""
    Input constants
"""
logMc=2.0                                       #log10 of chirp mass (Msun)
q=1.0                                           #mass ratio
#m1=36                                           #mass one (Msun)
#m2=30                                           #mass two (Msun)
#theta=0                                         #inclination angle (rad)
flow=10                                         #lower cutoff freq (Hz)
D=100       #standardisation, adjust distance in SNR calculator
#D=440                                           #distance (Mpc)
#default: m1=50, m2=50, D=0.778
#GW150914: m1=36, m2=30, D=440

#Constants derived from these used by both parts
#M = m1 + m2                                     #total mass
#eta = (m1*m2)/M**2                              #symmetric mass ratio
Mc=10**logMc                                      #chirp mass (Msun)
eta = q / ((1+q)**2)                            #substituting q in above eq.
M = Mc / (eta**(3/5))                           #standard chirp eq.

#reminder that we cannot calculate mergers with white dwarfs
if M*q/(1+q) <= 1.2:                            #white dwarf mass threshold
    raise ValueError('There is a white dwarf in this merger.')

"""
    Inspiral part
"""

#basic constants
pi=np.pi
c=2.99792458e8
#geometric unit conversion
Msunkg=1.9891e30
Msuns=4.923e-6
Msunkm=1.476

Dkm = D * 3.086e19                              #D in km

startx=(6.673e-11 * M*Msunkg*pi*flow*c**-3)**(2/3)              #eq 22
if M*q/(1+q) >= 3.0:                            #BH-BH mergers
    endx=(1/3)*(1 + (7/18)*eta)                     #eq 23, with 1/6 -> 1/3 (!)
else:                                           #BH-NS or NS-NS
    endx=(1/6)*(1 + (7/18)*eta)                     #eq 23, with 1/6
    #1/6 because 3RSch this is the ISCO for BH-NS and approx. the touching
    #point for NS-NS, so the end condition for both

def PNderiv(x):                                 #the diff eq for x
    #coefficients from Huerta article
    a0 = 12.8*eta
    a2 = eta*(-16*(924*eta + 743))/420
    a3 = 51.2*pi*eta
    a4 = (1903104*eta**2 + 3934368*eta + 1091296)*(eta/45360)
    #Buskirk also includes a8 to a12 (6PN)
    a8 = 170.799 - 742.551*eta + 370.173*eta**2 - 43.4703*eta**3 - \
        0.0249486*eta**4 + (14.143 - 150.692*eta)*np.log(x)
    a9 = 1047.25 - 2280.56*eta + 923.756*eta**2 + 22.7462*eta**3 - \
        102.446*np.log(x)
    a10 = 714.739 - 1936.48*eta + 3058.95*eta**2 - 514.288*eta**3 + \
        29.5523*eta**4 - 0.185941*eta**5 + (-3.00846 + 1019.71*eta + \
        1146.13*eta**2)*np.log(x)
    a11 = 3622.99 - 11498.7*eta + 12973.5*eta**2 - 1623*eta**3 + \
        25.5499*eta**4 + (83.1435 - 1893.65*eta)*np.log(x)
    a12 = 11583.1 - 45878.3*eta + 33371.8*eta**2 - 7650.04*eta**3 + \
        648.748*eta**4 - 14.5589*eta**5 - 0.0925075*eta**6 + (-1155.61 + \
        7001.79*eta - 2135.6*eta**2 - 2411.92*eta**3)*np.log(x) + \
        33.2307*np.log(x**2)
    Mdxdt = a0*x**5 + a2*x**6 + a3*x**6.5 + a4*x**7 + ((1/M)*(64/5)*eta*x**5 \
        * (1 + a8*x**4 + a9*x**4.5 + a10*x**5 + a11*x**5.5 + a12*x**6))
    return Mdxdt

dt=[]                                           #holds the timestep
xtimes=[0]                                      #times when x is evaluated
i=0                                             #loop counter
#Euler integration with variable timestep
x=[startx]                                      #post-Newtonian parameter
while x[i] <= endx:
    dt.append((10**-6*x[i])/PNderiv(x[i]))      #timestep with 10^-6 threshold
    x.append(x[i] + 10**-6*x[i])                #Euler step
    xtimes.append(xtimes[i] + dt[i])            #increment time-storer
    i += 1                                      #increment counter

i_time=np.zeros((len(xtimes)))                  #time in SI units
omega=np.zeros((len(xtimes)))                   #angular frequency
i_phase=np.zeros((len(xtimes)))                 #angular phase
freq=np.zeros((len(xtimes)))                    #frequency
r=np.zeros((len(xtimes)))                       #orbital radius
rdot=np.zeros((len(xtimes)))                    #analytical derivative of r
A1=np.zeros((len(xtimes)))                      #part of strain equations
A2=np.zeros((len(xtimes)))                      #part of strain equations
Aorth=np.zeros((len(xtimes)))                   #orthogonal polarisation
Adiag=np.zeros((len(xtimes)))                   #diagonal polarisation

for i in range(len(xtimes)):                    #real times
    i_time[i] = xtimes[i]*M*Msuns

#initial values
omega[0]=x[0]**1.5
freq[0]=omega[0]/(M*Msuns)

#for the phase integration, we are tied to the times at which we evaluated x
for i in range(len(xtimes)-1):                  #-1 because always i+1 terms
    omega[i+1] = x[i+1]**1.5
    i_phase[i+1] = i_phase[i] + omega[i+1]*dt[i]
    freq[i+1] = omega[i+1] / (M*Msuns*pi)
    
#plt.figure(1)
#plt.plot(i_time,phase)
    
#PN correction coefficients from B/BH article
r0pn = 1
r1pn = -1 + 0.333333*eta
r2pn = 4.75*eta + 0.111111*eta**2
r3pn = -7.51822*eta - 3.08333*eta**2 + 0.0246914*eta**3
for i in range(len(xtimes)):
    r[i] = r0pn*(1/x[i]) + r1pn + r2pn*x[i] + r3pn*x[i]**2
    rdot[i] = PNderiv(x[i]) * (-2*r0pn*x[i]**-2 + r2pn + 2*r3pn*x[i])
    
#plt.figure(1)
#plt.plot(i_time,r)

for i in range(len(xtimes)):                    #based on B/BH eq 9
    A1[i] = (-2*M*eta*(1/Dkm))*(rdot[i]**2 + (r[i]*omega[i])**2 + 1/r[i])
    A2[i] = (-2*M*eta*(1/Dkm))*(2*r[i]*rdot[i]*omega[i])

for i in range(len(xtimes)):
    Aorth[i] = A1[i]*np.cos(2*i_phase[i]) + A2[i]*np.sin(2*i_phase[i])
    Adiag[i] = A1[i]*np.sin(2*i_phase[i]) - A2[i]*np.cos(2*i_phase[i])

"""
    Merger part, up to frequency evolution
"""

if M*q/(1+q) >= 3.0:                            #merger part for BH-BH only
    Rscaled = D * 3.086e22
    
    sfin = 2*np.sqrt(3)*eta - (390/79)*eta**2 + (2379/287)*eta**3 - \
        (4621/276)*eta**4                           #final spin (eq 20)
    wqnm = 1 - 0.63*(1 - sfin)**0.3                 #quasi-normal modes (eq 19)
    
    #gIRS coefficients (Appendix C)
    Q = 2/((1 - sfin)**0.45)
    alpha = Q**-2 * (16313/562 + (21345/124)*eta)
    b = 16014/979 - (29132/1343)*eta**2
    C = 206/903 + (180/1141)*np.sqrt(eta) + (424/1205)*eta**2*(1/np.log(eta))
    kappa = 713/1056 - (23/193)*eta
    
    time = np.zeros((201))      #can't set it directly as range because of typing
    for i in range(201):
        time[i] = range(-100,101)[i]
        
    fhat = np.zeros((201))
    w = np.zeros((201))                             #really an omega (ang freq)
    for i in range(201):
        fhat[i] = (C/2) * (1 + 1/kappa)**(1 + kappa) * (1 - (1 + \
            (1/kappa)*np.exp(-2*time[i]*(1/b)))**-kappa)        #eq 17
        w[i] = 0.5 * wqnm * (1 - fhat[i])                       #eq 18 (note 0.5)
        
    fhatdot = np.zeros((201))
    #unusual averages at ends
    fhatdot[0] = fhat[1] - fhat[0]
    fhatdot[200] = fhat[200] - fhat[199]
    for i in range(1,200):                          #derivative as differences
        fhatdot[i] = 0.5*(fhat[i+1] - fhat[i-1])
        
    m_time = np.zeros((201))                        #time in s
    for i in range(201):
        m_time[i] = time[i]*M*4.923e-6

"""
    Inspiral array reduction
"""

#Attempt to fix bug with low-mass mergers
#(no longer necessary after MQ rewrite, but still useful for file size)

reduction_factor = 100                          #what factor to reduce by

#creating new inspiral arrays with every nth point
#only i_time and omega are used in the matching that was causing the bugs, but
#the others are used again later and all the inspiral arrays need to be the
#same length
reduced_i_time = [i_time[0]]
reduced_omega = [omega[0]]
reduced_freq = [freq[0]]                        #used by NS-NS instead of omega
reduced_i_phase = [i_phase[0]]
reduced_Aorth = [Aorth[0]]
reduced_Adiag = [Adiag[0]]
for i in range(reduction_factor,len(i_time),reduction_factor): #every nth point
    reduced_i_time.append(i_time[i])
    reduced_omega.append(omega[i])
    reduced_freq.append(freq[i])
    reduced_i_phase.append(i_phase[i])
    reduced_Aorth.append(Aorth[i])
    reduced_Adiag.append(Adiag[i])
    
#overwriting arrays with reduced versions
i_time = reduced_i_time
omega = reduced_omega
freq = reduced_freq
i_phase = reduced_i_phase
Aorth = reduced_Aorth
Adiag = reduced_Adiag
    
"""
    Matching frequency waveforms
"""

#The time-frequency pairs are (i_time,omega) inspiral (m_time,w) merger

if M*q/(1+q) >= 3.0:                            #merger part for BH-BH only
    for i in range(len(m_time)):
        m_time[i] += 0.04923      #having the merger start at t = 0 for convenience
        
    minMQ = 1000                                    #deliberately overly high
    MQ = []                                         #holds MQ values
    min_switch_ind = []
    
    #with this MQ method, we iterate for each frequency value in the merger part,
    #look for the closest frequency in the inspiral part, measure the df difference
    #between these points in the two waveforms, and then select the frequency with
    #the minimum df difference as the switching point, adjusting the waveform times
    #appropriately so the combined waveform is continuous in f and t
    
    def MQdiff(i):                                  #only one index needed here
        try:
            closest_index = np.searchsorted(omega, w[i], side='right')
            #with this method, we use searchsorted in the frequency domain instead
            #of the time domain
            #it does assume frequency increases monotonously, which it should
            df_diff = abs((w[i] - w[i-1])/(m_time[i] - m_time[i-1]) - \
                          (omega[closest_index] - omega[closest_index - 1])/ \
                          (i_time[closest_index] - i_time[closest_index - 1]))
        except:
            df_diff = np.nan
        #again, to get rid of potential IndexErrors for irrelevant frequencies
        
        return df_diff
    
    for i in range(1,len(m_time)):
        #MQ calculation at each point in merger waveform
        #EXCEPT 0, because then the [i-1] in MQdiff rolls over to the last value
        #in the array (and the first index shouldn't be the switch point anyway)
        MQ = MQdiff(i)
        if MQ < minMQ:
            minMQ = MQ                              #update minimum
            min_switch_ind = i                      #time of new minimum
    
    final_i_index = np.searchsorted(omega, w[min_switch_ind], side='right')
    #I guess MQdiff could be configured to return this somehow but doing it
    #separately is easier
    
    time_offset = i_time[final_i_index] - m_time[min_switch_ind]
    
    min_offset_m_time = np.zeros((len(m_time)))
    for i in range(len(m_time)):                    #offsetting to match i_time
        min_offset_m_time[i] = m_time[i] + time_offset
        
    #now we stitch the inspiral and merger frequency waveforms together    
    i_m_omega = []                                  #combined omega
    i_m_time = []                                   #combined time
    for i in range(final_i_index):                  #inspiral segment
        i_m_omega.append(omega[i])
        i_m_time.append(i_time[i])
        
    for i in range(min_switch_ind,len(m_time)):     #merger segment
        i_m_omega.append(w[i])
        i_m_time.append(min_offset_m_time[i])

"""
    Merger part phase and amplitude evolution
"""

#merger part phase calculation
m_phase = np.zeros((201 - min_switch_ind))      #orbital phase
#note: starts at min_switch_ind instead of the 0 used earlier for merger part
m_phase[0] = i_phase[final_i_index]             #matching phase of parts
for i in range(min_switch_ind + 1,201):
    m_phase[i - min_switch_ind] = m_phase[i - min_switch_ind - 1] + w[i]
                                                #Euler integration of eq 21
                                                
#stitching together phase
i_m_phase = np.zeros((final_i_index))           #combined phase
for i in range(final_i_index):                  #inspiral segment
    i_m_phase[i] = i_phase[i]
    
i_m_phase = np.concatenate((i_m_phase,m_phase)) #concatenating merger segment

#plt.figure(1)
#plt.plot(i_m_time,i_m_phase)

#merger part amplitude calculation
matching_amplitude = np.sqrt(Aorth[final_i_index]**2 + Adiag[final_i_index]**2)
#strain amplitude at end of inspiral segment
m_amp = np.zeros((201 - min_switch_ind))
for i in range(min_switch_ind,201):             #initial unscaled calculation
       m_amp[i - min_switch_ind] = (1/(2e0*w[i])) * ((abs(fhatdot[i]))/(1 + \
            alpha*(fhat[i]**2 - fhat[i]**4)))**(1/2) #* np.cos(m_phase[i - \
            #min_switch_ind])                    #eq 16 #with phase
            #(note 2e0)

scaling_ratio = matching_amplitude / m_amp[0]
for i in range(len(m_amp)):                     #rescaling for continuity
    m_amp[i] = scaling_ratio * m_amp[i]
    
#stitching together strain amplitude
i_m_amp = np.zeros((final_i_index))             #combined strain amplitude
#note: this does not include/differentiate the polarisations
for i in range(final_i_index):                  #inspiral segment
    i_m_amp[i] = np.sqrt(Aorth[i]**2 + Adiag[i]**2)
    
i_m_amp = np.concatenate((i_m_amp,m_amp))       #concatenating merger segment

#frequency in real units, for plotting purposes
i_m_freq = np.zeros((len(i_m_omega)))           #frequency in real units
for i in range(len(i_m_omega)):
    i_m_freq[i] = i_m_omega[i] / (M*Msuns*pi)   #omega -> freq conversion

"""
plt.figure(1)
plt.plot(i_m_time,i_m_amp)
plt.xlabel('Time (s)')
plt.ylabel('Strain amplitude')

plt.figure(2)
plt.plot(i_m_freq,i_m_amp)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Strain amplitude')
"""

"""
    Exporting array for SNR script to use
"""

#columns for exportarray: time, freq, amp

if M*q/(1+q) >= 3.0:                            #output for BH-BH
    exportarray = np.zeros((len(i_m_omega),3))
    for i in range(len(i_m_omega)):
        exportarray[i,0] = i_m_time[i]
        exportarray[i,1] = i_m_freq[i]
        exportarray[i,2] = i_m_amp[i]
        
else:                                           #output for BH-NS and NS-NS
    amp = np.zeros((len(Aorth)))                #replacement of i_m_amp
    for i in range(len(Aorth)):
        amp[i] = np.sqrt(Aorth[i]**2 + Adiag[i]**2)
    
    exportarray = np.zeros((len(i_time),3))    
    for i in range(len(i_time)):
        exportarray[i,0] = i_time[i]
        exportarray[i,1] = freq[i]
        exportarray[i,2] = amp[i]
    
#saving exportarray to a file
np.savetxt("output_" + str(logMc) + "_" + str(q) + ".dat",exportarray,\
           delimiter='\t',newline='\n')

"""    
dfarray = np.zeros((exportarray.shape[0]))
dfarray[0] = np.nan
for i in range(1,len(dfarray)):
    dfarray[i] = (exportarray[i,1] - exportarray[i-1,1]) / (exportarray[i,0] \
           - exportarray[i-1,0])
    #technically this is skewed to the left
 
plt.figure(1) 
plt.plot(exportarray[:,1],dfarray[:])
plt.xlabel('f (Hz)')
plt.ylabel('df (Hz/s)')
"""
