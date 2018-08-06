from imports import *
from truncate_cmap import *


def loadpickle(fname):
    fObj = open(fname, 'rb')
    self = pickle.load(fObj)
    fObj.close()
    return self


class TESSLC:

    def __init__(self, fname):
	try:
	    os.mkdir('TESSResults')
	except OSError:
	    pass
	self.fname = fname.replace('.fits','')
        self.fname_full = 'TESSResults/%s'%self.fname
        try:
            os.mkdir(self.fname_full)
        except OSError:
            pass
        self.pickleobject()


    def pickleobject(self):
        fObj = open('%s/TESSResults'%self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close() 
