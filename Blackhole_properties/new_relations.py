import numpy as np

# $\log (\frac{{M}_{BH}}{{M}_\odot})=2.85 \log (\frac{\sigma_0}{189}) + 1.16 \log (\frac{B/T}{0.437}) - 0.33 \log (\frac{\rho_{soi}}{601}) + 8.20$
def sigma_bt_rho_normalized_relation(log_sigma0,log_B_T,logRho_soi):
    a=2.85425535
    x1norm=2.2771162
    b=1.16151575
    x2norm=-0.35953965
    c=-0.3337952
    x3norm=2.77884623
    d=8.20260696
    return a*(log_sigma0)+b*(log_B_T)+c*(logRho_soi)+d-a*x1norm-b*x2norm-c*x3norm

#
def sigma_bt_rho_relation(log_sigma0,log_B_T,logRho_soi):
    a=2.77679962
    b=1.16427157
    c=-0.36642744
    d=3.32060273
    return a*log_sigma0+b*log_B_T+c*logRho_soi+d

# $\log (\frac{{M}_{BH}}{{M}_\odot})=2.57 \log \sigma_0 + 0.38 (\log R_e - \log \rho_{soi}) + bvc + 2.61$
def sigma_r_rho_bvc_relation(log_sigma0,log_R_e_sph_eq_kpc,logRho_soi_approx_new,bvtc):
    a=2.5665281441619349
    b=0.37680597349565176
    c=-0.37680597349565176
    d=1.0
    e=2.6084185089278386
    return a*(log_sigma0)+b*(log_R_e_sph_eq_kpc)+c*(logRho_soi_approx_new)+d*(bvtc)+e

# $\log (\frac{{M}_{BH}}{{M}_\odot})=2.57 \log \sigma_0 + \log(B/T) - 0.43 \log \rho_{soi} - 0.24 \ {Pseudobulge} + 4.01$
def sigma_bt_rho_pseudobulge_relation(log_sigma0,log_BT,logRho_soi_approx_new,pseudobulge):
    a=2.5665281441619349
    b=1.0
    c=-0.4279452531210189
    d=-0.2393228544761872
    e=4.0131431262964635
    return a*(log_sigma0)+b*(log_BT)+c*(logRho_soi_approx_new)+d*(pseudobulge)+e

# $\log (\frac{{M}_{BH}}{{M}_\odot})=3 \log \sigma_0 - (\log \rho_{soi})^{0.58} + e^{\log(B/T)-{Pseudobulge}} + 2.55$
def sigma_rho_bt_pseudobulge_relation(sigma,rho,bt,pseudobulge):
    a=3.0
    b=0.575925501030858
    c=2.547821423461232
    return a*sigma-rho**b+np.exp(bt-pseudobulge)+c

# $\log (\frac{{M}_{BH}}{{M}_\odot})=0.91 \log {M}^{*}_{{sph}} - 0.41(\log \rho_{soi}+{Pseudobulge}) - 0.27$
def Msph_rho_pseudobulge_relation(Msph,rho,pseudobulge):
    a=0.9074312354176397
    b=-0.406124683495805
    c=-0.2718555890272426
    return a*Msph+b*(rho+pseudobulge)+c

'''# RMS too high
def sigma_Msph_BT_relation(sigma,Msph,BT):
    a=0.58404844826419349
    return sigma+a*(Msph+BT)
'''

# $\log (\frac{{M}_{BH}}{{M}_\odot})= \log L_{B}-\frac{1.73}{bvc}$
def blum_bvc_relation(blum,bvc):
    a=1.727784796655826
    return blum-(a/bvc)

# $\log (\frac{{M}_{BH}}{{M}_\odot})= 0.44^{bvc} + \log L_{B} - 5.52$
def bvc_blum_relation(bvc,blum):
    a=4.444375083380111
    b=-5.52239593792091
    return a**bvc+blum+b

# $\log (\frac{{M}_{BH}}{{M}_\odot})=\log \sigma_0 (3.61 + 0.23({Core}-{Pseudobulge}))$
def sigma_core_pseudobulge(sigma,core,pseudobulge):
    a=0.23472318279549
    b=3.612014128255
    return sigma*(a*(core-pseudobulge)+b)

# $\log (\frac{{M}_{BH}}{{M}_\odot})=4.14 \log \sigma_0 + 0.45 ({Core}-{Pseudobulge})-1.14$
def sigma_core_pseudobulge2(sigma,core,pseudobulge):
    a=4.135926771730385
    b=0.45114601637524224
    c=-1.1422301553820672
    return a*sigma+b*(core-pseudobulge)+c

# $\log (\frac{{M}_{BH}}{{M}_\odot})= \log {M}^{*}_{{sph}} + 0.36 ({Core}-{Pseudobulge}) -2.42$
def Msph_core_pseudobulge(Msph,core,pseudobulge):
    a=0.36492546046497276
    b=-2.419965528706549
    return Msph+a*(core-pseudobulge)+b

# $\log (\frac{{M}_{BH}}{{M}_\odot})= \log {M}^{*}_{{sph}} + 0.33 ({Core}-{Pseudobulge}) + 0.17 bvc -2.54$
def Msph_core_pseudobulge_bvc(Msph,core,pseudobulge,bvc):
    a=0.3346068919534675
    b=0.17447573368509578
    c=-2.5420715074368143
    return Msph+a*(core-pseudobulge)+b*bvc+c

# $\log (\frac{{M}_{BH}}{{M}_\odot})= \log R_e + 0.67 {Core} + bvc + 7.06$
def Re_core_bvc_relation(Re,core,bvc):
    a=0.66638191376642007
    b=7.063105852055458
    return Re+a*core+bvc+b


# spiral/LTG galaxy
# $\log (\frac{{M}_{BH}}{{M}_\odot})= (\log \sigma_0 + 1.40)(\logv_{max}-\tan\phi)$
def sigma_phi_v_relation(log_sigma0,phi,v):
    a=1.3969453234456826
    return (log_sigma0+a)*(v-phi)

# spiral/LTG galaxy
# $\log (\frac{{M}_{BH}}{{M}_\odot})= 3.89 (\log \sigma_0 - \tan\phi)$
def sigma_phi_relation(sigma,phi):
    a=3.8894438058884226
    b=-3.8894438058884226
    return a*sigma+b*phi

# spiral/LTG galaxy
# $\log (\frac{{M}_{BH}}{{M}_\odot})= 1.51 \log v_{max}(\log v_{max}-\tan\phi) $
def phi_v_relation(phi,v):
    a=1.514256269494091
    return a*v*(v-phi)

# spiral/LTG galaxy
# $\log (\frac{{M}_{BH}}{{M}_\odot}) = \frac{8.93+\frac{0.61 \tan\phi}{M^*_{gal}(M^*_{gal}-\frac{2.50}{\tan\phi})}}{\tan\phi+0.98}$
def mgal_phi_relation(mgal,phi):
    a=8.928946755726638
    b=0.6121230044044399
    c=2.5005639649031903
    d=0.9796146336168521
    return (a+(b*phi/(mgal*(mgal-(c/phi)))))/(phi+d)

# Easy-to-use relation
# $\log (\frac{{M}_{BH}}{{M}_\odot}) = \log \sigma_0 + \log \textup{M}^{*}_{\textup{sph}} - 0.42 \ \textup{Pseudobulge} - 4.59$
def sigma_Msph_pseudobulge_relation(sigma,Msph,pseudobulge):
    a=-0.4187824174909826
    b=-4.5895796522317545
    return sigma+Msph+a*pseudobulge+b

# Easy-to-use relation
# $\log (\frac{{M}_{BH}}{{M}_\odot}) = 0.93 (\log \sigma_0)^2 + 0.56 \log R_{e,sph,maj} + 3.2)$
def sigma_Re_maj_relation(sigma,Re):
    a=0.9271154736686101
    b=0.5607269304350105
    c=3.2029708446894385
    return a*sigma**2+b*Re+c

# Easy-to-use relation
# $\log (\frac{{M}_{BH}}{{M}_\odot}) = 3.59 \log \sigma_0 + 0.50 \log R_{e,sph,maj} - 0.50 \textup{Pseudobulge}$
def sigma_Re_maj_Pseudobulge_relation(sigma,Re,pseudobulge):
    a=3.589249391739994
    b=0.50273659760562763
    c=-0.50273659760562763
    return a*sigma+b*Re+c*pseudobulge