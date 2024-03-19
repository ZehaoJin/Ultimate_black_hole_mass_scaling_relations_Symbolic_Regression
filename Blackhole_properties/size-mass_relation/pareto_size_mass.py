# produce a dense pareto front of black hole scaling relations
import numpy as np
import pandas as pd
from pysr import PySRRegressor
import re
from tqdm import tqdm

# params
# already did low_adv
para_set = 'low'  # low/easy
operator_set = 'sim'  #adv/sim

df_full = pd.read_csv('SMBH_Data_03_07_24.csv',header=1)

# paras = ['M*_sph','M*_gal','R_e_sph_maj','log_R_e_sph_maj','R_e_sph_eq_kpc','log_R_e_sph_eq_kpc',
#         'R10', 'logR10', 'logR10phi', 'Rh', 'logRh', 'logRhphi', 'logRemajphi', 'logReeqphi']

paras = ['log_sigma0','M*_gal','logR10','logRh']

if paras[-1]!='M_BH':
    paras.append('M_BH')
paras.append('M_BH_std_sym')

obs=df_full.copy()
obs = obs[paras].dropna(axis='index',how='any')
print(len(obs))

y = obs['M_BH'].to_numpy()
w = 1/obs['M_BH_std_sym'].to_numpy()**2

X = obs.iloc[:,:-2].to_numpy()

# start over 50 times, each time evolve 500 generations
#evolutions = 50
evolutions = 5
niterations = 5000


denoise=False
ncyclesperiteration=2000
adaptive_parsimony_scaling=20
verbosity=0
return_model=True
maxsize=25
optimizer_algorithm="NelderMead"
optimizer_iterations=100
optimizer_nrestarts=5
should_simplify=True
temp_equation_file=True  # turn this on will prevent the auto-save of hall_of_fame
tempdir='/data/zj448/SR/Ultimate_paper/temp/'  # dummy dir to save temp files and then being auto deleted

for epoch in tqdm(range(evolutions),desc='evolution'):

    warm_start=True
    model = PySRRegressor(
        binary_operators=["+", "-", "*", "/"],
        warm_start=warm_start,
        denoise=denoise,
        niterations=1,
        ncyclesperiteration=ncyclesperiteration,
        adaptive_parsimony_scaling=adaptive_parsimony_scaling,
        verbosity=verbosity,
        precision=64,
        maxsize=maxsize,
        optimizer_algorithm=optimizer_algorithm,
        optimizer_iterations=optimizer_iterations,
        optimizer_nrestarts=optimizer_nrestarts,
        should_simplify=should_simplify,
        temp_equation_file=temp_equation_file,
        tempdir=tempdir,
        procs=64
        )

    for iter in range(niterations):
        

        model.fit(X=X, y=y, weights=w)


        equations = model.equations_
    
        number_matching_pattern = r"(?<![a-zA-Z0-9_.])[+-]?(\d+\.\d+|\.\d+|\d+\.|\d+)(?:[eE][-+]?\d+)?"

        variable_matching_pattern = r'\bx\d+'


        # Count number of constants:
        equations["number_constants"] = [len(re.findall(number_matching_pattern, eq)) for eq in equations["equation"]]

        # count number of (unique) variables
        equations["variables"]=[set(re.findall(variable_matching_pattern, eq)) for eq in equations["equation"]]
        equations["number_variables"] = [len(re.findall(variable_matching_pattern, eq)) for eq in equations["equation"]]
        equations["unique_number_variables"] = [len(set(re.findall(variable_matching_pattern, eq))) for eq in equations["equation"]]
        

        # Compute log likelihood (for example)
        #equations["log_like"] = - equations["loss"] * len(X)
        #equations["log_like"] = -0.5 * len(X) * np.log(2*np.pi) - (np.log(obs['M_BH_std_sym'])).sum() - 0.5 * equations["loss"] * w.sum()

        # Compute AIC:
        #equations["aic"] = 2 * equations["number_constants"] - 2 * equations["log_like"]

        equations['evolutions']=epoch
        equations['iterations']=iter

        if epoch==0 and iter==0:
            equations.to_csv('/data/zj448/SR/Ultimate_paper/pareto_archive/pareto_sm.csv',index=False)
        else:
            equations.to_csv('/data/zj448/SR/Ultimate_paper/pareto_archive/pareto_sm.csv',index=False,mode='a',header=False)
