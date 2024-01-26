# produce a dense pareto front of black hole scaling relations
import numpy as np
import pandas as pd
from pysr import PySRRegressor
import re
from tqdm import tqdm

low_scatter_para=['ETG','T-type','Bar', 'Disk', 'Ring', 'Core', 'Multiple', 'Compactness', 'AGN',
       'Pseudobulge', 'BCG', 'cD','M*_sph', 'M*_gal', 'log_B/T',
       'log_sigma0', 'log_R_e_sph_maj','log_R_e_sph_eq_kpc', 'log_n_sph_maj', 'log_n_sph_eq', 'log(I_e,sph,maj/M_Sun/pc^2)',
       'log(I_e,sph,eq/M_Sun/pc^2)', 'Concentration_Index',
       'avg_Rho_1kpc_Exact_All', 'r1_density_approx', 'log10(R10_kpc)',
       'logRho_R10_approx', 'log_rho10_Exact', 'log10(R90_kpc)',
       'logRho_R90_approx', 'log_rho_90_Exact_all', 'Rho_re_spatial',
       'SR_pc_All', 'Rho_SR_pc_All', 'CR_def1_approx_new',
       'Rho_cr_def1_approx_new', 'CR_def2_approx_new',
       'Rho_CR_def2_approx_new', 'Sr(pc)_2_using_Falserm_drho',
       'Log_Approx_Avg_density_10pc', 'log_Rho_e_Exact_new',
       'logRho_e_approx_New', 'logRho_soi_approx_new',
       'log_Rho_soi_exact_new', 'Avg_Rho_Re_Exact_all',
       'Avg_Rho_soi_exact_all', 'Avg_Rho_re_Exact_all', 'Rho_re_Exact_all',
       'Rho_r_soi_2BH_approx', 'Log_Avg_Rho_10kpc_approx',
       'Log_Avg_Rho_10kpc_exact_final', 'Log_Avg_Rho_100pc_approx',
       'Log_Avg_Rho_5kpc_approx', 'Log_Avg_rho_5kpc_exact_all', 'ube', 'bve',
       'dc', 'bvtc', 'bri25', 'mabs', 'blum', 'logblum', 'logSigma0sph',
       'LogSigma0', 'R10', 'logR10', 'logR10phi', 'Rh', 'logRh', 'logRhphi',
       'logHalo','B-V','V-[3.6]','GJC23W1-W2','GJC23W2-W3','GJC23log(M*,gal/M_sun)',
       'GJC23log(SFR)','GJC23log(sSFR)','M_BH']

easy_obs_para=['LogSigma0','Concentration_Index','logSigma0sph','log_sigma0','dc','logRhphi','M*_sph','ube','bri25','bve','bvtc','logR10phi','M*_gal','log_B/T',
 'logRh','log_n_sph_eq','blum','log_R_e_sph_maj','logblum','log_n_sph_maj','logR10','Pseudobulge','AGN','Multiple','Ring','BCG','Disk','cD',
 'Bar','Core','Compactness','ETG','T-type','M_BH','log10(R10_kpc)','log10(R90_kpc)','B-V','V-[3.6]','GJC23W1-W2','GJC23W2-W3','GJC23log(M*,gal/M_sun)','M_BH']


df_full = pd.read_csv('SMBH_Data_01_22_24.csv',header=1)


paras=low_scatter_para.copy()
#paras=easy_obs_para.copy()

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
        binary_operators=["+", "-", "*", "/","pow"],
        unary_operators=["exp","log10"],
        constraints={"pow": (9, 1)}, # power laws can have 9 complexity left argument, but only 1 complexity in the right argument.
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
        equations["log_like"] = - equations["loss"] * len(X)

        # Compute AIC:
        equations["aic"] = 2 * equations["number_constants"] - 2 * equations["log_like"]

        equations['evolutions']=epoch
        equations['iterations']=iter

        if epoch==0 and iter==0:
            equations.to_csv('/data/zj448/SR/Ultimate_paper/pareto_archive/pareto.csv',index=False)
        else:
            equations.to_csv('/data/zj448/SR/Ultimate_paper/pareto_archive/pareto.csv',index=False,mode='a',header=False)
