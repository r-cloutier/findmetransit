import numpy as np
import os

start, num = 0, 0
NP, Nrps = 15, 6
Ps = np.logspace(np.log10(.5), np.log10(351), NP)
#rpRss = np.linspace(.1, .3, NrpRs)
rps = np.logspace(np.log10(.5), np.log10(15), Nrps)

for i in range(NP):
    for j in range(Nrps):
        f = open('jobscript_sensv4', 'r')
        g = f.read()
        f.close()

        g = g.replace('<<fnum>>', '%i'%(start+num))
        num += 1
        g = g.replace('<<P>>', '%.5f'%Ps[i])
        g = g.replace('<<rp>>', '%.5f'%rps[j])
        h = open('jobscript', 'w')
        h.write(g)
        h.close()

        #os.system('cat jobscript')
        os.system('qsub jobscript')
        os.system('rm jobscript')
