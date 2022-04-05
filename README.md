# curemodels
Penalized mixture cure models for high-dimensional data

``EM.R`` and ``gmifs.R`` implement the Weibull MCM estimation methods developed by Han et al. Competing methods are implemented in ``Weibull.R``, ``cmix.R`` and ``glmnet_cox.R``. The data generating processes are implemented in ``data_generation.R`` and evaluation metrics are included in ``evaluation.R``. 

``simul_all_in_one.R`` is a main function where the other functions are called and different methods are compared with each other.
