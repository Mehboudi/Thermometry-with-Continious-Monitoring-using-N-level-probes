import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad, odeint, solve_ivp
from scipy.linalg import logm, expm

col_palette = ["navy","slateblue","blueviolet","purple","deeppink","firebrick","coral","peru","darkorange","goldenrod","olive","darkolivegreen","seagreen","mediumturquoise"]

plt.rcParams.update({
    "text.latex.preamble": r"\usepackage{amsmath}",
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8

})

# Parameters
s = 1.         #   ohmicity
Omega = 10.0     #   cut-off
gamma = 0.1     #   strength
beta  = 1.0     #   temperature

w0  = 1.0       #   spin-frequency

# Spectral density
def J(w,s):
    return np.pi/2 * gamma * w**s / Omega**(1-s) * np.exp(-w/Omega)

# Bosonic occupation number
def n(w,beta):
    return np.exp(-beta*w)/(1-np.exp(-beta*w))

# True decay rates
def Gamma10(w0,t,beta,s):
    def integrand(w):
        return J(w,s) / (2 * np.pi) * t * (np.sinc((w+w0)*t/(np.pi)) * (n(w,beta) + 1) + np.sinc((w-w0)*t / (np.pi)) * n(w,beta))
    return quad(integrand,0,1)[0] + quad(integrand,1,np.infty,limit=1000)[0]

def Gamma01(w0,t,beta,s):
    def integrand(w):
        return J(w,s) / (2 * np.pi) * t * (np.sinc((w+w0)*t/(np.pi)) * (n(w,beta)) + np.sinc((w-w0)*t / (np.pi)) * (n(w,beta)+1))
    return quad(integrand,0,1)[0] + quad(integrand,1,np.infty,limit=1000)[0]

def ME(t,p,w0,beta,s):
    G10 = Gamma10(w0,t,beta,s)
    G01 = Gamma01(w0,t,beta,s)
    Gamma = np.array([[-G01,G10],
                      [G01,-G10]])
    return Gamma.dot(p)

def numericalRate(w0,tau,beta,s,n_eval=200):
    p00 = np.array([1.0,0.0])
    p01 = np.array([0.0,1.0])

    res0 = solve_ivp(ME,[0.0,tau],p00,args=(w0,beta,s),method="LSODA",t_eval=np.linspace(0,tau,n_eval))
    res1 = solve_ivp(ME,[0.0,tau],p01,args=(w0,beta,s),method="LSODA",t_eval=np.linspace(0,tau,n_eval))

    pt0 = res0.y
    t0  = res0.t
    pt1 = res1.y
    t1  = res1.t

    M = np.transpose(np.vstack((pt0,pt1)).reshape((2,2,len(t0))),(1,0,2))
    # M = np.vstack((pt0,pt1)).reshape((2,2,len(t0)))

    return pt0, t0, pt1, t1, M

def obtainTrueFishy(w0,tau,beta,s,n_eval=200,dBeta=1e-3):
    TMP = numericalRate(w0,tau,beta,s,n_eval=n_eval)
    times = TMP[1]
    M0 = TMP[-1][:,:,1:]
    M1 = numericalRate(w0,tau,beta+dBeta,s,n_eval=n_eval)[-1][:,:,1:]
    M_1 = numericalRate(w0,tau,beta-dBeta,s,n_eval=n_eval)[-1][:,:,1:]

    p0 = M0[0,1,:] / (1 - M0[0,1,:] + M0[0,0,:])
    p1 = 1 - p0

    if np.sum(p0<0) + np.sum(p1<0) > 0:
        print("FUCK")

    Fishy = np.zeros((n_eval))

    dM = (M1 - M0)/(dBeta)

    Fishy = 0.0
    for j in range(2):
        Fishy += p0 * dM[j,0,:]**2 / M0[j,0,:] + p1 * dM[j,1,:]**2 / M0[j,1,:]

    return times, np.concatenate((np.array([0.0]),Fishy / times[1:]))

def obtainThermalFishy(w0,beta):
    Z = 1 + np.exp(-beta*w0)
    p0 = 1/Z
    p1 = np.exp(-beta*w0)/Z

    avgE = p1*w0

    return p0 * (0 - avgE)**2 + p1 * (w0 - avgE)**2

# Effective decay rates
def Gamma10_eff(w0,tau,beta,s):
    def integrand(w):
        return J(w,s) * tau/(4*np.pi) * (np.sinc((w+w0)*tau/(2*np.pi))**2 * (n(w,beta) + 1) + np.sinc((w-w0)*tau / (2*np.pi))**2 * n(w,beta))
    return quad(integrand,0,1,limit=1000)[0] + quad(integrand,1,np.infty,limit=1000)[0]

def Gamma01_eff(w0,tau,beta,s):
    def integrand(w):
        return J(w,s) * tau/(4*np.pi) * (np.sinc((w-w0)*tau/(2*np.pi))**2 * (n(w,beta) + 1) + np.sinc((w+w0)*tau / (2*np.pi))**2 * n(w,beta))
    return quad(integrand,0,1,limit=1000)[0] + quad(integrand,1,np.infty,limit=1000)[0]

# Beta derivatives of effective decay rates
def DBeta_Gamma10_eff(w0,tau,beta,s):
    def integrand(w):
        return - w * J(w,s) * 1/(4*np.pi) * tau * (np.sinc((w-w0)*tau/(2*np.pi))**2 + np.sinc((w+w0)*tau / (2*np.pi))**2) * (np.exp(-beta*w)/(1-np.exp(-beta*w))**2)
    return quad(integrand,0,1,limit=1000)[0] + quad(integrand,1,np.infty,limit=1000)[0]

def DBeta_Gamma01_eff(w0,tau,beta,s):
    def integrand(w):
        return - w * J(w,s) * 1/(4*np.pi) * tau * (np.sinc((w+w0)*tau/(2*np.pi))**2 + np.sinc((w-w0)*tau / (2*np.pi))**2) * (np.exp(-beta*w)/(1-np.exp(-beta*w))**2)
    return quad(integrand,0,1,limit=1000)[0] + quad(integrand,1,np.infty,limit=1000)[0]

# Obtain Fisher information for the effective rate model
def obtainEffFishyMarkov(w0,tau,beta,s,n_eval=200):
    times = np.linspace(1e-3,tau,n_eval)
    g10_range = np.vectorize(lambda t : Gamma10_eff(w0,t,beta,s))(times)
    g01_range = np.vectorize(lambda t : Gamma01_eff(w0,t,beta,s))(times)
    DBeta_g10_range = np.vectorize(lambda t : DBeta_Gamma10_eff(w0,t,beta,s))(times)
    DBeta_g01_range = DBeta_g10_range # np.vectorize(lambda t : DBeta_Gamma01_eff(w0,t,beta,s))(times)

    # Obtain stationary state
    p1 = g10_range / (g01_range + g10_range)
    p0 = 1 - p1

    # Obtain Fisher
    return times, (p0 * (DBeta_g10_range**2 / g10_range) + p1 * (DBeta_g01_range**2 / g01_range))

# Obtain Fisher information for the effective rate model
def obtainEffFishy(w0,tau,beta,s,n_eval=200,dBeta=1e-3):
    times = np.linspace(1e-3,tau,n_eval)

    # Obtain M(beta)
    g10_range = np.vectorize(lambda t : Gamma10_eff(w0,t,beta,s))(times)
    g01_range = np.vectorize(lambda t : Gamma01_eff(w0,t,beta,s))(times)

    Geff = (np.array([[-g01_range,g10_range],
                      [g01_range,-g10_range]])*times).transpose(2,0,1)
    
    M0 = expm(Geff).transpose(1,2,0)

    p0 = M0[0,1,:] / (1 - M0[0,1,:] + M0[0,0,:])
    p1 = 1 - p0

    # Obtain M(beta+dBeta)
    g10_range = np.vectorize(lambda t : Gamma10_eff(w0,t,beta+dBeta,s))(times)
    g01_range = np.vectorize(lambda t : Gamma01_eff(w0,t,beta+dBeta,s))(times)

    Geff = (np.array([[-g01_range,g10_range],
                      [g01_range,-g10_range]])*times).transpose(2,0,1)
    
    M1 = expm(Geff).transpose(1,2,0)

    Fishy = np.zeros((n_eval))

    dM = (M1 - M0)/(dBeta)

    Fishy = 0.0
    for j in range(2):
        Fishy += p0 * dM[j,0,:]**2 / M0[j,0,:] + p1 * dM[j,1,:]**2 / M0[j,1,:]

    return times, Fishy / times

###################
# KIDS PLAYGROUND #
###################
# tau = 10
# n_eval = 10
# t_range = np.linspace(0.0,tau,n_eval)

# M = numericalRate(w0,tau,beta,s,n_eval=n_eval)[-1]

# g10_range = np.vectorize(lambda t : Gamma10_eff(w0,t,beta,s))(t_range)
# g01_range = np.vectorize(lambda t : Gamma01_eff(w0,t,beta,s))(t_range)

# plt.figure()
# plt.plot(t_range,g01_range*t_range,label="10")
# plt.plot(t_range,M[1,0,:],label="01")
# plt.legend(loc="best")
# plt.show()

#########################
# RATIO versus TIME BIN #
#########################
t_range = np.linspace(1e-3,20,100)
if False:
    plt.figure(figsize=(6,4),dpi=300)
    for k,s in enumerate([0.5]):
        g10_range = np.vectorize(lambda t : Gamma10_eff(w0,t,beta,s))(t_range)
        g01_range = np.vectorize(lambda t : Gamma01_eff(w0,t,beta,s))(t_range)
        plt.plot(t_range,g10_range/g01_range,color=col_palette[2*k],label=r"$\Gamma_{10}^{\rm eff} / \Gamma_{01}^{\rm eff}$, $s="+str(s)+r"$")

        g10_true, g01_true = np.vectorize(lambda t : numericalRate(w0,t,beta,s))(t_range)
        plt.plot(t_range,g10_true/g01_true,color=col_palette[2*k+1],label=r"$\Gamma_{10}^{\rm true} / \Gamma_{01}^{\rm true}$, $s="+str(s)+r"$")
    plt.hlines(np.exp(beta*w0),min(t_range),max(t_range),color="black",linestyle="--",label=r"detailed balance $e^{\beta\omega_0}$")
    plt.xlabel(r"time bin width $\Delta t = [\omega_0^{-1}]$")
    plt.ylabel(r"ratio of rates $\Gamma_{\rm out} / \Gamma_{\rm in} = [1]$")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("figures/ratio_vs_time.jpg",dpi=300)
    plt.show()
    plt.close()

###############################
# FISHER INFO versus TIME BIN #
###############################
if True:
    tau=20
    n_eval=100
    t_range = np.linspace(0.0,tau,n_eval)

    thermalFishy = obtainThermalFishy(w0,beta)

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(7,2.5))

    ax1.hlines(thermalFishy,min(t_range),max(t_range),color="black",linestyle=":",label=r"$F_\beta^\mathrm{thermal}$")

    lines=["--","-",":"]
    for k,s in enumerate([0.75,1.0]):
        times_true, FishyTrue = obtainTrueFishy(w0,tau,beta,s,n_eval)
        # times_effMarkov, FishyEffMarkov = obtainEffFishyMarkov(w0,tau,beta,s,n_eval)
        times_eff, FishyEff = obtainEffFishy(w0,tau,beta,s,n_eval)
        ax1.plot(times_true,FishyTrue*times_true,linestyle=lines[k],color=col_palette[2*k],label=r"$F_\beta^\mathrm{true}$, $s="+str(s)+r"$")
        ax1.plot(times_eff,FishyEff*times_eff,linestyle=lines[k],color=col_palette[2*k+1],label=r"$F_\beta^\mathrm{eff}$, $s="+str(s)+r"$")

        ax2.plot(times_true,FishyTrue,linestyle=lines[k],color=col_palette[2*k],label=r"$F_\beta^\mathrm{true}$, $s="+str(s)+r"$")
        # ax2.plot(times_effMarkov,FishyEffMarkov,linestyle=lines[k],color=col_palette[2*k],label=r"$F_\beta^\mathrm{eff,markov}$, $s="+str(s)+r"$")
        ax2.plot(times_eff,FishyEff,linestyle=lines[k],color=col_palette[2*k+1],label=r"$F_\beta^\mathrm{eff}$, $s="+str(s)+r"$")

    ax1.set_xlabel(r"time bin $\Delta = [1/\omega_0]$")
    ax1.set_ylabel(r"FI per step")
    ax1.legend(loc="best")

    ax2.set_xlabel(r"time bin $\Delta = [1/\omega_0]$")
    ax2.set_ylabel(r"FI rate")
    ax2.legend(loc="best")
    plt.tight_layout()
    plt.savefig("figures/fishy_fisher.jpg",dpi=300)
    plt.show()
    plt.close()

#####################
# RATIO versus BETA #
#####################
beta_range = np.linspace(1e-3,5,20)
if False:
    plt.figure(figsize=(6,4))
    for k,tau in enumerate([0.1,1.0,10.0,100.0]):
        g10_range = np.vectorize(lambda b : Gamma10_eff(w0,tau,b,0.5))(beta_range)
        g01_range = np.vectorize(lambda b : Gamma01_eff(w0,tau,b,0.5))(beta_range)
        plt.semilogy(beta_range,g10_range/g01_range,color=col_palette[k],linestyle="--",marker=".",label=r"$\Gamma_{10}^{\rm eff} / \Gamma_{01}^{\rm eff}$, $s=0.5$, $\Delta \tau="+str(tau)+r"/\omega_0$")
    plt.plot(beta_range,np.exp(beta_range*w0),color="black",linestyle="--",label=r"detailed balance $e^{\beta\omega_0}$")
    plt.xlabel(r"inverse temperature $\beta = [\omega_0^{-1}]$")
    plt.ylabel(r"ratio of rates $\Gamma_{\rm out} / \Gamma_{\rm in} = [1]$")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("figures/ratio_vs_beta.jpg",dpi=300)
    plt.close()
    # plt.show()

####################
# Populations plot #
####################
if True:
    # t_range = np.linspace(0,10,10)
    fig, (ax1,ax2) = plt.subplots(2,1)
    t_range = np.linspace(1e-3,20,300)
    ax1.plot(t_range,np.vectorize(lambda t : Gamma10(w0,t,beta,1.0))(t_range),label="true, s=1")
    ax1.plot(t_range,np.vectorize(lambda t : Gamma10(w0,t,beta,0.5))(t_range),label="true, s=0.5")
    ax1.plot(t_range,np.vectorize(lambda t : Gamma10_eff(w0,t,beta,1.0))(t_range),label="eff, s=1")
    ax1.plot(t_range,np.vectorize(lambda t : Gamma10_eff(w0,t,beta,0.5))(t_range),label="eff, s=0.5")
    ax1.legend(loc="best")
    # plt.legend(loc="best")
    # plt.plot(t_range,np.vectorize(lambda t : Gamma01(w0,t,beta,2.0))(t_range),label="s=2")
    pt0, t0, pt1, t1, M = numericalRate(w0,max(t_range),beta,.5)
    ax2.plot(t0,pt0[0,:],label="pt0_0")
    # plt.plot(t0,pt0[1,:],label="pt0_1")
    # pt0, t0, pt1, t1 = numericalRate(w0,50,beta,1.0)
    # plt.plot(t0,pt0[0,:],label="pt0_0")
    # plt.plot(t0,pt0[1,:],label="pt0_1")
    # plt.plot(t1,pt1[0,:],label="pt1_0")
    ax2.plot(t1,pt1[1,:],label="pt1_1")
    ax2.legend(loc="best")
    plt.show()

