def m_sigma_relation(log_sigma0):
    a=5.82
    b=-5.221994574764372
    return a*(log_sigma0)+b


def sigma_phi_v_relation(log_sigma0,phi,v):
    a=1.3969453234456826
    return (log_sigma0+a)*(v-phi)


def sigma_phi_relation(sigma,phi):
    a=3.8894438058884226
    b=-3.8894438058884226
    return a*sigma+b*phi


def phi_v_relation(phi,v):
    a=1.514256269494091
    return a*v*(v-phi)


def mgal_phi_relation(mgal,phi):
    a=8.928946755726638
    b=0.6121230044044399
    c=2.5005639649031903
    d=0.9796146336168521
    return (a+(b*phi/(mgal*(mgal-(c/phi)))))/(phi+d)
