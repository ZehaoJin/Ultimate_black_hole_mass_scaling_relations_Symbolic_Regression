low_scatter_para=['ETG', 'T-type', 'Bar', 'Disk', 'Ring', 'Core', 'Multiple',
    'Compactness', 'AGN', 'Pseudobulge', 'BCG', 'cD', 
    'M*_sph', 'M*_gal', 'log_B/T', 'log_sigma0', 'log_R_e_sph_maj', 'log_R_e_sph_eq_kpc',
    'log_n_sph_maj', 'log(I_e,sph,maj/M_Sun/pc^2)',
    'log10(R10_kpc)', 'logRho_R10_approx', 'log_rho10_Exact', 'log10(R90_kpc)', 'logRho_R90_approx',
    'log_rho_90_Exact_all', 'Log_Approx_Avg_density_10pc', 'log_Rho_e_Exact_new', 'logRho_e_approx_New',
    'logRho_soi_approx_new', 'log_Rho_soi_exact_new', 'Log_Avg_Rho_10kpc_approx', 
    'Log_Avg_Rho_10kpc_exact_final', 'Log_Avg_Rho_100pc_approx', 'Log_Avg_Rho_5kpc_approx',
    'Log_Avg_rho_5kpc_exact_all', 'dc', 'bvtc', 'mabs', 'blum', 'logblum', 'logSigma0sph',
    'LogSigma0', 'R10', 'logR10', 'logR10phi', 'Rh', 'logRh', 'logRhphi', 'logHalo', 'GJC23W1-W2',
    'GJC23W2-W3', 'GJC23log(M*,gal/M_sun)', 'GJC23log(SFR)', 'GJC23log(sSFR)', 'log<Sigma>_e',
    'log<Sigma>_h', 'M_BH']

easy_obs_para=['LogSigma0', 'logSigma0sph', 'log_sigma0', 'dc', 'logRhphi', 'M*_sph', 'bvtc',
    'logR10phi', 'M*_gal', 'log_B/T', 'logRh', 'blum', 'log_R_e_sph_maj',
    'logblum', 'log_n_sph_maj', 'logR10', 'Pseudobulge', 'AGN', 'Multiple', 'Ring', 'BCG',
    'Disk', 'cD', 'Bar', 'Core', 'Compactness', 'ETG', 'T-type', 'log10(R10_kpc)', 'log10(R90_kpc)',
    'GJC23W1-W2', 'GJC23W2-W3', 'GJC23log(M*,gal/M_sun)', 'M_BH']

low_scatter_para_std=low_scatter_para.copy()
for i in low_scatter_para:
    if i!='M_BH':
        low_scatter_para_std.append(i+'_std')
low_scatter_para_std.append('M_BH_std_sym')

easy_obs_para_std=easy_obs_para.copy()
for i in easy_obs_para:
    if i!='M_BH':
        easy_obs_para_std.append(i+'_std')
easy_obs_para_std.append('M_BH_std_sym')