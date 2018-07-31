    params, resultsGPfin, mufin, sigfin = joint_LC_fit(params, thetaGPout, bjd, f, fcorr, ef)

def joint_LC_fit(params, thetaGP, bjd, f, fcorr, ef):
    '''Iteratively optimize the planetary model parameters
