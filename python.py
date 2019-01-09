# Title: My Fifth KMC Code
# Date: 1/9/2018
# Author: Jian Ren Lim

from scipy.integrate import odeint    # import odeint for calculating analytical solution
import random,math,os,matplotlib.pyplot as plt,numpy as np   # import random number generator, math, and plot functions
path = 'JR KMC output txt_file'     # create a folder to contain the output txt.files
if not os.path.exists(path):        # check if this directory exist. If it does, we do not need to create this directory
    os.makedirs(path)               # create directory

constants = [1.8,3.2]       # rate constants iteration for  k1, and k2
num = [100,500,2000]       # number of initial A molecules
numB = [0, 200, 1000]       # number of initial B molecules
for k1 in constants:  # k1 will either be 1.8 or 3.2
    for k2 in constants:  # k2 will either be 1.8 or 3.2
        for N in num:  # N (number of initial A molecules) will either be 100, 1000, or 4000
            for NB in numB:     # NB (number of initial B molecules) will either be 0, 500, 2000

                # Section A: Analytical Part

                p0 = [N, NB, 0]  # number of initial A, B, and C molecules. C will always be zero to start with.


                def model(p, t):  # definition of 'model' function for analytical solution
                    dxdt = -k1 * p[0]
                    dydt = k1 * p[0] - k2 * p[1]
                    dzdt = k2 * p[1]
                    dpdt = [dxdt, dydt, dzdt]  # collect dA/dt, dB/dt, and dC/dt into placeholder dp/dt
                    return dpdt  # return dp/dt


                t = np.linspace(0, 5)  # let time span from 0 to 5 (e.g. seconds)
                p = odeint(model, p0, t)  # get analytical solution for 0<t<5

                # Section B: KMC part

                T = []     # declare empty array T for time
                A = []     # declare empty array A for # of a
                B = []     # declare empty array B for # of b
                C = []     # declare empty array C for # of c
                time = 0    # let t0=0 to start with
                a = p0[0]   # declare initial number of a exactly the same as section A
                b = p0[1]   # declare initial number of b exactly the same as section A
                c = p0[2]   # declare initial number of c exactly the same as section A

                # Create txt.file to store KMC data
                text_file = open('%s/k1=%.1f_k2=%.1f,NA0=%d,NB0=%d' % (path,k1, k2, N, NB),"w")     # create a new txt.file
                text_file.write("%-14s%-8s%-8s%-8s%-8s%-8s" % ('t', 'proc', 'A', 'B', 'C', 'rand')) # prepare labels for the data we're going to input later
                # 'proc' shows the process taking place while 'rand' shows the random number between 0 and rtot that determines which rxn will take place

                while c < N+NB:  # let time loop until all molecules have been depleted except C
                    r1 = k1 * a  # rate expression of A -> B
                    r2 = k2 * b  # rate expression of B -> C
                    rtot = r1 + r2  # add rate 1 & rate 2 for the purpose of estimating dt (i.e. how fast the next reaction will occur)

                    rand = random.random() * (rtot)  # select random number 'rand' between 0 and rtot
                    dtime = -math.log(
                        random.random()) / rtot  # let the multiple of timestep and instantenous rtot be the negative logarithm of a random number: r*t=-log(n) or n=exp(-r*t)
                    T.append(time)      # append time into T array
                    A.append(a)         # append # of a into A array
                    B.append(b)         # append # of b into B array
                    C.append(c)         # append # of c into C array
                    time = time + dtime  # update time for the following loop

                    # When selecting which reaction should take place, we look at the proportionality of each rate to the total rate.
                    # The higher the proportionality, the higher the chance that this reaction will take place in this timestep.
                    # However, due to random collision of molecules, we will assume a random number to be the barrier for reaction to take place.
                    # If the reaction is smaller/larger than one of the rates, then that reaction will take place, and vice versa.
                    # The arragement (e.g. rand < r2, rand >r2, rand <r1, rand >r1) doesn't matter as long as it is consistent throughout the code.

                    if (rand < r2):  # if r2 is greater than rand, r2 takes place
                        b = b - 1;  # b is being used by 1 unit
                        c = c + 1  # c is being produced by 1 unit
                        text_file.write("\n%.8f    b->c    %-8s%-8s%-8s%.2f" % (time, a, b, c, rand))   # input data after this reaction

                    else:  # if r2 is smaller than rand, r2 takes place
                        a = a - 1;  # a is being used by 1 unit
                        b = b + 1  # b is being produced by 1 unit
                        text_file.write("\n%.8f    a->b    %-8s%-8s%-8s%.2f" % (time, a, b, c, rand))   # input data after this reaction

                plt.plot(t, p[:, 0], 'b--', label=r'$Conc.A_{anal}$')  # this function plot A from analytical part
                plt.plot(t, p[:, 1], 'g--', label=r'$Conc.B_{anal}$')  # plot B
                plt.plot(t, p[:, 2], 'r--', label=r'$Conc.C_{anal}$')  # plot C
                plt.plot(T, A, label='$A_{kmc}$')  # this function plot A from kmc part
                plt.plot(T, B, label='$B_{kmc}$')  # plot B
                plt.plot(T, C, label='$C_{kmc}$')  # plot C
                plt.title('$k_1= %.1f s^{-1},k_2= %.1f s^{-1},N_{A0}= %d,N_{B0}= %d$' % (k1, k2, N, NB));  # Show the title of individual plot
                plt.ylabel('# of molecules')  # label on y-axis
                plt.xlabel('time')  # label on x-axis
                plt.legend(loc='best')  # show the legend at location that doesn't obstruct line plots
                plt.show()  # show the plot we've made

                text_file.close()   # close the current txt.file

