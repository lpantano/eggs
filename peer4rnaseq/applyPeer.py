
import peer
import scipy as SP
import pylab as PL
import pdb


print "Simple PEER application. All default prior values are set explicitly as demonstration."
y = SP.loadtxt("/home/lpantano/scripts/reproducibility/peer4rnaseq/tsiexp.txt",delimiter=",")
K = 8
Nmax_iterations = 100
model = peer.PEER()

# set data and parameters
model.setNk(K) #number of factor for learning
model.setPhenoMean(y) # data for inference
# set priors (these are the default settings of PEER)
model.setPriorAlpha(0.001,0.1);
model.setPriorEps(0.1,10.);
model.setNmax_iterations(Nmax_iterations)
# perform inference
model.update()

#covs = SP.loadtxt("data/covs.csv", delimiter=",")
#model.setCovariates(covs2) # covariates (e.g batch, RNA quality, ...)
#model.update()

