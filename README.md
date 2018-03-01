# Dark matter induced asymmetry in atomic clock data, and annual modulation

Benjamin M. Roberts

2018-02-05, 16:23

Companion code to paper:

_Precision measurement noise asymmetry and its annual modulation as a dark matter signature_,
B. M. Roberts & A. Derevianko


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/benroberts999/DM-ClockAssymetry/blob/master/LICENSE)


## simulateDMclock.py [python3]

Program simulates (white noise) atomic clock data [s(1)],
and injects DM signals at random times, according to Poisson distribution.

The assumptions from the paper are not assumed in the simulation.
For example, the crossing is allowed to occur at any time, not just in between
data sample times.

After generating random white noise for the s(1) clock data, we inject DM events
at random times following a Poisson distribution with average rate, R_0.
We calculate each signal chi for Gaussian monopoles by drawing v and rho
randomly according to their respective probability distributions (see paper).
Then, the skew is calculated from the data, without any knowledge of the
original noise distribution, or the assumed DM parameters,
confirming that it does indeed scale as stated in the paper, and also confirms
the statement that the equation slightly underestimates the skew.

Calculates the total skew, and compares to approximate theoretical value.
Theory underestimates the skew by factor of ~3, justifying the 'conservative'
approximation.

Will output the over-all skew information to the screen, does not write to file.

Will also calculate the annual modulation in the skew (broken up by week),
so long as 'num_years' greater than 1.
The skew for each week is writen to file, in the form:
  *  week    k3(week)     dk3     N_M

It then perform the FFT, and returns the calculated modulation amplitude for
the skew.
Note: Just reutns for the MAX of the FFT - so not necisarily the correct
frequency!


## numericallyCalculateSkew [c++]

Calculates the expected skew numerically, for a range of x0 and r0t0 values.
Uses Gaussian monopole, and full standard halo model distributions.

Values are saved to a text file, in form: x0 r0t0 skew
Default filename: numericalSkew.txt
Note: will just over-ride the file if it already exists.

p_s(s) = r0t0 * Integrate[ pc(r) * px(s-r) , {r, -infty, infty} ]
        + (1-r0t0) * pc(s),

where pc is the 'intrinsic' clock noise (Gaussian).
px (p_chi) is the DM signal distribution (in absence of noise). See paper.
