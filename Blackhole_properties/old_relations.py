import numpy as np
'''my fit
def m_sigma_relation(log_sigma0):
    a=5.39254382
    b=-4.0642592
    return a*(log_sigma0)+b
'''

# https://arxiv.org/pdf/1908.06838.pdf
# general
# $\log (\frac{{M}_{BH}}{{M}_\odot})= 6.10 \log (\frac{\sigma_0}{200})+8.27$
def m_sigma_relation(log_sigma0):
    a=6.10
    b=8.27
    return a*(log_sigma0-np.log10(200.))+b

'''didn't work
# https://arxiv.org/pdf/1810.04888.pdf
# general
# $\log (\frac{{M}_{BH}}{{M}_\odot})=2.41 \log [\frac{\log(B/T)}{-0.77}] + 7.15$
def m_bt_relation(log_bt):
    a=2.41
    b=7.15
    c=-0.77
    return a*np.log10(log_bt/c)+b
'''

# Ultimate_black_hole_mass_scaling_relations_Symbolic_Regression/Blackhole_properties/reference_relations_0502.ipynb
# general
# $\log (\frac{{M}_{BH}}{{M}_\odot})=2.40 \log(B/T) + 9.11$
def m_bt_relation(log_bt):
    a=2.40365436
    b=9.11049
    return a*log_bt+b


# Ultimate_black_hole_mass_scaling_relations_Symbolic_Regression/Blackhole_properties/reference_relations_0502.ipynb
# general
# $\log (\frac{{M}_{BH}}{{M}_\odot})=-0.93 \log \rho_{soi} + 10.64$
def m_rho_soi_relation(log_rho_soi):
    a=-0.93004258
    b=10.63954454
    return a*log_rho_soi+b

# Ultimate_black_hole_mass_scaling_relations_Symbolic_Regression/Blackhole_properties/reference_relations_0502.ipynb
# general
# $\log (\frac{{M}_{BH}}{{M}_\odot})=1.31 \log M^*_{sph} - 5.83$
def m_msph_relation(log_m_sph):
    a=1.31446409
    b=-5.83301935
    return a*log_m_sph+b

# Ultimate_black_hole_mass_scaling_relations_Symbolic_Regression/Blackhole_properties/reference_relations_0502.ipynb
# general
# $\log (\frac{{M}_{BH}}{{M}_\odot})=1.94 \log L_{B} -12.07603801$
def m_blum_relation(logblum):
    a=1.94161099
    b=-12.07603801
    return a*logblum+b

# Ultimate_black_hole_mass_scaling_relations_Symbolic_Regression/Blackhole_properties/reference_relations_0502.ipynb
# general
# $\log (\frac{{M}_{BH}}{{M}_\odot})=1.63 \log R_{e,sph,eq} + 7.68$
def m_Re_eq_relation(log_Re):
    a=1.63073002
    b=7.67555795
    return a*log_Re+b

# Ultimate_black_hole_mass_scaling_relations_Symbolic_Regression/Blackhole_properties/reference_relations_0502.ipynb
# general
# $\log (\frac{{M}_{BH}}{{M}_\odot})=1.51 \log R_{e,sph,maj} + 7.61$
def m_Re_maj_relation(log_Re):
    a=1.50980505
    b=7.60613021
    return a*log_Re+b


# https://arxiv.org/pdf/1908.06838.pdf
# LTG galaxies
# $\log (\frac{{M}_{BH}}{{M}_\odot})= 5.82 \log (\frac{\sigma_0}{200})+8.17$
def m_sigma_IMBH_relation(log_sigma0):
    a=5.82
    b=8.17
    return a*(log_sigma0-np.log10(200.))+b

# https://arxiv.org/pdf/1707.04001.pdf
# LTG galaxies
# $\log (\frac{{M}_{BH}}{{M}_\odot})= 7.01-0.171(\phi-15^\circ)$
def m_phi_relation(tan_phi):
    a=7.01
    b=0.171
    c=15
    phi=np.rad2deg(np.arctan(tan_phi))
    return a-b*(phi-c)

# https://arxiv.org/pdf/1901.06509.pdf
# sprial galaxies
# $\log (\frac{{M}_{BH}}{{M}_\odot})= 10.62 \log (\frac{v_{max}}{210})+7.22$
def m_vmax_relation(log_vmax):
    a=10.62
    b=7.22
    return a*(log_vmax-np.log10(210.))+b

# https://arxiv.org/pdf/2209.14526.pdf
# sprial galaxies
# $\log (\frac{{M}_{BH}}{{M}_\odot})=3.23 \log (M^*_{sph}/10^{11}) + 7.91$
def m_mgal_relation(log_mgal):
    a=3.23
    b=7.91
    return a*(log_mgal-np.log10(1e11))+b
