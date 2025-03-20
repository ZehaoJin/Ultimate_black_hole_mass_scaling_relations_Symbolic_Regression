from general_fit_function import general_fit

import pandas as pd
from global_parameter import *
from tqdm import tqdm
import numpy as np


## parameters
## low/easy: low_scatter/easy_obs, adv/sim: operators
# already did low_adv, low_sim, easy_adv, easy_sim
filename = 'pareto_easy_sim_refit.csv'  
t_eq=pd.read_csv('/data/zj448/SR/Ultimate_paper/pareto_archive/'+filename)

df_full = pd.read_csv('SMBH_Data_03_06_24.csv',header=1)

if 'low' in filename:
    paras = low_scatter_para
    paras_std = low_scatter_para_std
if 'easy' in filename:
    paras = easy_obs_para
    paras_std = easy_obs_para_std

booleans=['ETG','Bar','Disk','Ring','Core','Multiple','Compactness','AGN','Pseudobulge','BCG','cD']
for b in booleans:
    df_full[b+'_std']=0
    
df=df_full[paras_std].dropna(axis='index',how='any').copy()


for i, column in enumerate(paras):
    if column != 'M_BH':
        df.rename(columns={column: 'x'+str(i), column+'_std': 'x'+str(i)+'_std'}, inplace=True)
    else:
        df.rename(columns={'M_BH': 'y', 'M_BH_std_sym': 'y_std'}, inplace=True)


t_eq['final_equation']=pd.Series(dtype='object')
t_eq['logdet']=0.
t_eq['sqrt_cov']=pd.Series(dtype='object')
t_eq['final_loss']=0.
t_eq['optimized_params']=pd.Series(dtype='object')
t_eq['N_para']=0.
t_eq['residual_array']=pd.Series(dtype='object')
t_eq['ad_p_value_array']=pd.Series(dtype='object')
t_eq['final_BIC']=0.
t_eq['final_BIC_i']=0.



for row in tqdm(t_eq.iterrows(),total=len(t_eq)):
    if np.isnan(row[1]['LLL']):
        equation = row[1]['equation']
        intrinsic_scatter_guess = None
        simplification = True
    else:
        equation = row[1]['refit_equation']
        intrinsic_scatter_guess = row[1]['intrinsic_scatter']
        simplification = False

    try:
        result = general_fit(equation, df, initial_lr=0.05,
                            verbose=False, simplification=simplification,
                            intrinsic_scatter_guess=intrinsic_scatter_guess)
    except ValueError:
        try:
            result = general_fit(equation, df, initial_lr=0.01,
                                verbose=False, simplification=simplification,
                                intrinsic_scatter_guess=intrinsic_scatter_guess)
        except ValueError:
            result = {'fitted_equation': np.nan, 'logdet': np.nan, 'sqrt_cov': np.nan, 'loss': np.nan, 'optimized_params': np.nan,
                    'N_para': np.nan, 'residual_array': np.nan, 'ad_p_value_array': np.nan, 'BIC': np.nan, 'BIC_i': np.nan}

        
    t_eq.at[row[0],'final_equation'] = result['fitted_equation']
    t_eq.at[row[0],'logdet'] = result['logdet']
    t_eq.at[row[0],'sqrt_cov'] = result['sqrt_cov']
    t_eq.at[row[0],'final_loss'] = result['loss']
    t_eq.at[row[0],'optimized_params'] = result['optimized_params']
    t_eq.at[row[0],'N_para'] = result['N_para']
    t_eq.at[row[0],'residual_array'] = result['residual_array']
    t_eq.at[row[0],'ad_p_value_array'] = result['ad_p_value_array']
    t_eq.at[row[0],'final_BIC'] = result['BIC']
    t_eq.at[row[0],'final_BIC_i'] = result['BIC_i']


    t_eq.to_csv('/data/zj448/SR/Ultimate_paper/pareto_archive/'+filename[:-10]+'_generalfit.csv',index=False)