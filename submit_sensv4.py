import numpy as np
import os

start, num = 0, 0
NP, Nrps = 11, 8
Ndays_field = 27.4
Pgrid = np.logspace(np.log10(.5), np.log10(Ndays_field), NP+1)
rpgrid = np.logspace(np.log10(.5), np.log10(15), Nrps+1)

for i in range(NP):
    for j in range(Nrps):
        f = open('jobscript_sensv4', 'r')
        g = f.read()
        f.close()

        g = g.replace('<<fnum>>', '%i'%(start+num))
        num += 1
        g = g.replace('<<Pin>>', '%.5f'%Pgrid[i])
        g = g.replace('<<Pout>>', '%.5f'%Pgrid[i+1])
        g = g.replace('<<rpin>>', '%.5f'%rpgrid[j])
        g = g.replace('<<rpout>>', '%.5f'%rpgrid[j+1])
	g = g.replace('<<Ndays>>', '%.2f'%Ndays_field)
        h = open('jobscript', 'w')
        h.write(g)
        h.close()

        #os.system('cat jobscript')
        os.system('qsub jobscript')
        os.system('rm jobscript')
