import numpy as np
import time


A = np.random.rand(64,1000000)
toremove = [7,10,11,13,14,15,19,22,23,25,26,27,28,29,30,31,34,35,37,38,39,40,41,42,43,44,45,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
start=time.time();
np.delete(A,toremove,0);
end=time.time();
print "Time: ", (end-start)*1000.0
