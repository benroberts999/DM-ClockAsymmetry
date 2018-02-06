#!/usr/bin/python3
import random
import math
import numpy as np
from numpy import sqrt, sin, cos, pi
from scipy.special import erfinv, erf
from scipy import stats
from scipy.stats import binned_statistic
################################################################################
# Benjamin M. Roberts
# 2018-02-05, 16:23
#
# Program simulates (white noise) atomic clock data [s(1)],
# and injects DM signals at random times, according to Poisson distribution.
#
# The assumptions from the paper are not assumed in the simulation.
# For example, the crossing is allowed to occur at any time, not just in between
# data sample times.
#
# After generating random white noise for the s(1) clock data, we inject DM events
# at random times following a Poisson distribution with average rate, $R_0$.
# We calculate each signal chi for Gaussian monopoles by drawing v and rho
# randomly according to their respective probability distributions (see paper).
# Then, the skew is calculated from the data, without any knowledge of the
# original noise distribution, or the assumed DM parameters,
# confirming that it does indeed scale as stated in the paper, and also confirms
# the statement that the equation slightly underestimates the skew.
#
# Calculates the total skew, and compares to approximate theoretical value.
# Theory underestimates the skew by factor of ~3, justifying the 'conservative'
# approximation.
#
# Will output the over-all skew information to the screen, does not write to file.
#
# Will also calculate the annual modulation in the skew (broken up by week),
# so long as 'num_years' greater than 1.
# The skew for each week is writen to file, in the form
#   week    k3(week)     dk3     N_M
#
# It then perform the FFT, and returns the calculated modulation amplitude for
# the skew.
# Note: Just reutns for the MAX of the FFT - so not necisarily the correct
# frequency!
#
################################################################################
# Input parameters:

# number of 'years' to generate data for [nb: 0.02yrs =~ 1 week]
# If >1, will look for annual modulation.
# Annual modulation works best for integers!
num_years=0.02

# Data sample period (in seconds). ~1 is typical.
tau0 = 1.

# Standard deviation of intrinsic clock noise [always 1 - just sets 'units']
sd=1.

# DM event parameters:
# Typical DM signal magnitude (chi_0)
x0=1.0*sd
# expected average rate PER SAMPLE:
R0tau0=0.005
#DM object size/width (in km). Note: should be <<10**3
#(OK if larger, but skew will be spuriously small, due to hidden d factor in x0)
d=1.
# Velocity modulation \delta v / v. dv_on_v=0.05 according to SHM
dv_on_v=0.05  #shouldn't change this
################################################################################

# Global constants - don't chang:
DAYSINYEAR=365.25   # days in a year (average)
MINSINDAY =1440     # seconds in a day
SECSINDAY =86400     # seconds in a day

#-------------------------------------------------------------------------------
def generateClockNoise(num_points,sd):
    # Function that returns a numpy array of simulate clock (frequency) data
    # Assumes white noise. standard deviation and number of points input.
    print("Generating white frequency clock noise for",num_points,"sample points...")
    x0=0.   # example centre for Gaussian
    u_lst = [random.uniform(0,1) for _ in range (num_points)]
    clk_data = [(x0+sqrt(2.)*sd*erfinv(2*u_lst[i]-1)) for i in range (num_points)]
    print("...done")
    return np.array(clk_data)
#end generateClockNoise

#-------------------------------------------------------------------------------
def signal_s1(j,t0,tau0,v,d,rho):
    # DM monopole signal (Gaussian profile)
    # Integrated (S1) signal - for timing data (not freq.)
    # j: time 'now'. j0: impact time. tau0: sample time
    # v,d,rho: velocity, width, impact parameter
    # NOTE: does not include 'x0' (magnitude) [i.e. x0=1]
    # Chi_0 is taken out. Means, should return 1 if rho=0.7d, v=300.
    # NB: this means d has also been factored out! stored in chi_0!
    x = v/d
    arg1 = x*(j-t0)*tau0
    arg2 = x*(j-t0-1)*tau0
    sig = erf(arg1) - erf(arg2)
    impact_factor = math.exp(-math.pow(rho/d,2))
    fac = 300.*1.77245/(2*v)  # New: factor our Gamma_0
    return fac*sig*impact_factor
#end signal

#-------------------------------------------------------------------------------
def inject_DM(R0tau0,d,x0,dv_on_v):
#step through each epoch.
#Probability of an event occuring at that epoch depends on
#   - event rate
#   - phase (time of year)
    print("Injecting DM events: R0t0=",R0tau0,"x0=",x0,"d=",d,"dv/v=",dv_on_v," ...")
    num=0   # counts the number of injected events
    window = int(2*d/(100.*tau0)) + 2 # +/- this # of data points for each signal
    for j0 in range(num_points):
        phase = (j0%pts_yr)/pts_yr
        #"p" is average rate each SAMPLE PERIOD [expected # events in tau_0]
        p = R0tau0*(1+dv_on_v*cos(2*pi*phase))
        # Number of "events" during this sample period:
        events_epoch = np.random.poisson(p) # Use Poisson distribution!
        for n in range(events_epoch):
            # randomly choose t0 (anywhere within [j0-1,j0] epoch):
            t0 = j0 - random.uniform(0,1)
            #Randomly choose rho (impact parameter); in range [0,d], r^2
            rho = d*sqrt(random.uniform(0,1))
            # Randomly choose v [v. good Gaussian approx. to SHM vel. disto!]
            u_v = random.uniform(0.1,0.9)
            v0 = 360.*(1+dv_on_v*cos(2*pi*phase)) #check!! sd change also??
            v = v0 + 1.41*130.*erfinv(2*u_v-1)
            # work out signal for each 'epoch' in window
            for j in range (j0-window,j0+window+1):
                if j<0 or j>=num_points: continue
                clk_data[j]+=x0*signal_s1(j,t0,tau0,v,d,rho)
            num+=1 # count number of events
    #END inject events
    print("...done\n")
    return num
# END inject_DM()


################################################################################
################################################################################

print("\n ~~ Simulateing",num_years,"years =",num_years*52.179,"weeks ~~\n")

# Generate clock noise:
num_points = int(num_years*DAYSINYEAR*SECSINDAY/tau0)
clk_data = generateClockNoise(num_points,sd)
# Number of data points per year (needed for 'phase'):
pts_yr = int(num_points/num_years)

#Inject DM parameters
num = inject_DM(R0tau0,d,x0,dv_on_v)

print("Events injected/expected:",num,"/",int(R0tau0*num_points))
print("Num points per year:",pts_yr,"Number of years:",num_years)

# re-zero the data: (Not really needed)
old_mean = np.mean(clk_data)
clk_data -= old_mean
print("mean (before):",old_mean,", (after):",np.mean(clk_data))
print("std deviation:",np.std(clk_data),"\n")

# Calculate "skewness" of data:
print("Overal results for entire data set:")
k3=stats.skew(clk_data)
print("Skew of data:",k3)
dk3 = math.sqrt(6./num_points)
kth = (2/5)*R0tau0*math.pow(x0/sd,3) #small r0t0 expansion
print("k3 Theory=   ",kth)
print("Theory/Calc= ",kth/k3)
print("sqrt(6/n)=   ",dk3)
print("k3/dk3=      ",k3/dk3)
b=stats.skewtest(clk_data)
print(b)
print("")

# Look for annual modulations
# Only does if more than 1 year is simulated
if(num_years>1.):
    file = open("testfile.txt", "w")
    points_in_week = int(num_points/num_years/52.179)
    num_weeks = 1+int(num_points/points_in_week)
    k3week_list = []
    file.write("# week k3(week)  dk3   N_M")
    for iw in range(num_weeks):
        tj_i = points_in_week*iw
        tj_f = tj_i + points_in_week
        if (tj_f>num_points): tj_f = num_points
        N_M = tj_f-tj_i
        week_data = clk_data[tj_i:tj_f]
        k3_week = stats.skew(week_data)
        dk3_week = math.sqrt(6./N_M)
        k3week_list += [k3_week]
        outlst=[iw,k3_week,dk3_week,N_M]
        ostr = "\n"
        for el in outlst:
            ostr = ostr+str(el)+" "
        file.write(ostr)
    k3week_mean = np.mean(k3week_list)
    print("Weekly k3 mean=",k3week_mean)
    k3week_list -= k3week_mean
    abs_fft = np.absolute(np.fft.fft(k3week_list))
    #print(abs_fft)
    half_f = int(math.ceil((len(abs_fft)-1)/2.))
    #print(len(abs_fft),half_f)
    maxfft = np.amax(abs_fft[0:half_f])
    maxfreq = np.argmax(abs_fft[0:half_f])
    print("k3^(m)=",2.*maxfft/num_weeks,"(for the ",maxfreq,"cycles-per-T_total frequency)\n")
