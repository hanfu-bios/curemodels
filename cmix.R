
################## cmix ################## 

library(reticulate)
use_python("~/.conda/envs/cmix_env/bin/python", required=TRUE)

py_string="
import sys
sys.path.append('~/c-mix/C-mix/')
import numpy as np
import pylab as pl
import pandas as pd
from QNEM.inference import QNEM
from QNEM.simulation import CensoredGeomMixtureRegression
from sklearn.model_selection import ShuffleSplit
from lifelines.utils import concordance_index as c_index_score
from time import time
from matplotlib import rc
X=r.X
Y = np.array(r.Y)
delta = np.array(r.delta).astype(int)

features_names = range(X.shape[1])
n_samples, n_features = X.shape

tol = 1e-6
eta = 0
fit_intercept = True
gamma_chosen = '1se'
warm_start = True
grid_size = 30
metric = 'C-index'
verbose = True

model = 'CURE'  # 'C-mix', 'CURE'

learner = QNEM(l_elastic_net=0., eta=eta, max_iter=100, tol=tol,\
  warm_start=warm_start, verbose=verbose, model=model,\
  fit_intercept=fit_intercept)
learner.n_features = n_features

## Cross-validation ##
learner.cross_validate(X, Y, delta, n_folds=5, verbose=False, eta=eta,\
  grid_size=grid_size, metric=metric)
avg_scores = learner.scores.mean(axis=1)
l_elastic_net_best = learner.l_elastic_net_best
if gamma_chosen == '1se':
  l_elastic_net_chosen = learner.l_elastic_net_chosen
if gamma_chosen == 'min':
  l_elastic_net_chosen = l_elastic_net_best

grid_elastic_net = learner.grid_elastic_net  

## Run selected model with l_elasticNet_chosen ##
learner = QNEM(l_elastic_net=l_elastic_net_chosen, eta=eta, tol=tol,\
  warm_start=warm_start, verbose=verbose, model=model,\
  fit_intercept=fit_intercept)
learner.n_features = n_features
learner.fit(X, Y, delta)

coeffs = learner.coeffs
"
