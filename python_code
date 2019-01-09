# Title: My Third KMC Code
# Date: 1/3/2018
# Author: Jian Ren Lim

from scipy.integrate import odeint    # import odeint for calculating analytical solution
import random,math,matplotlib.pyplot as plt,numpy as np   # import random number generator, math, and plot functions

constants = [1,10]       # rate constants iteration for volume, k1, and k2
num = [100,1000,4000]    # number of initial A molecules
for k1 in constants:  # k1 will only be either 1 or 10
    for k2 in constants:  # k2 will only be either 1 or 10
        for N in num:  # N (number of initial A molecules) will either be 100, 1000, or 4000

            # Section A: Analytical Part

            p0 = [N, 0, 0]  # number of initial A, B, and C molecules


            def model(p, t):  # definition of 'model' function for analytical solution
                dxdt = -k1 * p[0]
                dydt = k1 * p[0] - k2 * p[1]
                dzdt = k2 * p[1]
                dpdt = [dxdt, dydt, dzdt]  # collect dA/dt, dB/dt, and dC/dt into placeholder dp/dt
                return dpdt  # return dp/dt


            t = np.linspace(0, 5)  # let time span from 0 to 5 (e.g. seconds)
            p = odeint(model, p0, t)  # get analytical solution for 0<t<5

            # Section B: KMC part

            T = [];
            A = [];
            B = [];
            C = [];
            time = 0  # declare empty array of T for time , # of A, B, C. Let t0=0
            a = p0[0];
            b = p0[1];
            c = p0[2]  # declare initial number of a, b, c exactly the same as section A

            while c < N:  # let time loop until all molecules have been depleted except C
                r1 = k1 * a;  # rate of A -> B
                r2 = k2 * b;  # rate of B -> C
                rtot = r1 + r2;  # sum rate 1 & 2 for the purpose of estimating how fast the next reaction occur

                rand = random.random() * (rtot);  # select random number 'rand' between 0 and rtot
                dtime = -math.log(
                    random.random()) / rtot  # multiple of timestep and instantenous rtot is equals to the negative logarithm of a random constant: n=exp(-r*t)
                T.append(time);
                A.append(a);
                B.append(b);
                C.append(c);  # append values into array
                time = time + dtime  # update time for the following loop

                # When selecting which reaction should take place, we look at the proportionality of each rate to the total rate.
                # The higher the proportionality, the higher the chance that this reaction will take place in this timestep.
                # However, due to random collision of molecules, we will assume a random number to be the barrier for reaction to take place.
                # If the reaction is smaller/larger than one of the rates, then that reaction will take place, and vice versa.
                # The arragement (e.g. rand < r2, rand >r2, rand <r1, rand >r1) doesn't matter as long as its consistent throughout.

                if (rand < r2):  # if r2 is greater than rand, r2 takes place
                    b = b - 1;  # b is being used by 1 unit
                    c = c + 1  # c is being produced by 1 unit
                else:  # if r2 is smaller than rand, r2 takes place
                    a = a - 1;  # a is being used by 1 unit
                    b = b + 1  # b is being produced by 1 unit

            plt.plot(t, p[:, 0], 'b--', label=r'$Conc.A_{anal}$')  # this function plot A from analytical part
            plt.plot(t, p[:, 1], 'g--', label=r'$Conc.B_{anal}$')  # plot B
            plt.plot(t, p[:, 2], 'r--', label=r'$Conc.C_{anal}$')  # plot C
            plt.plot(T, A, label='$A_{kmc}$')  # this function plot A from kmc part
            plt.plot(T, B, label='$B_{kmc}$')  # plot B
            plt.plot(T, C, label='$C_{kmc}$')  # plot C
            plt.title(
                '$k_1= %d s^{-1},k_2= %d s^{-1},N_{A0}= %d $' % (k1, k2, N));  # Show the title of individual plot
            plt.ylabel('# of molecules')  # label on y-axis
            plt.xlabel('time')  # label on x-axis
            plt.legend(loc='best')  # show the legend at location that doesn't obstruct line plots
            plt.show()  # show the plot we've made
