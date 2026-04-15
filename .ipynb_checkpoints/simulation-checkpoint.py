import numpy as np
from helper_functions import *
import sys

def main():
    # parameters
    L = 1
    Nsites = int(sys.argv[1])
    h = L / Nsites

    Dn = 0.8
    Dc = 1
    gamma = 2
    kappa = 1
    lambd = 25

    # homogeneous solution (unsteady)
    rhon = 2
    rhoc = rhon * gamma / kappa
    mu = Nsites / rhoc
    mu = 10

    Deff = (Dn - lambd * (T1(rhoc, mu, 70) - T0(rhoc, mu)))
    print(f'Deff: {Deff}')
    mu, rhoc

    Dn = Dn / h ** 2
    Dc = Dc / h ** 2
    lambd = lambd / h ** 2

    occupancies_n, dict_n = initialize_fields(Nsites, 1 / rhon, rho=rhon)
    occupancies_c, dict_c = initialize_fields(Nsites, 2 / rhoc, rho=rhoc)
    chemo_rates, cum_chemo_rates = initialize_chemo_rates(occupancies_n, occupancies_c, mu)

    N, Nc = int(np.sum(occupancies_n)), int(np.sum(occupancies_c))

    occupancies = np.vstack((occupancies_n, np.zeros_like(occupancies_n), occupancies_c))
    dicts = [dict_n, None, dict_c]

    T = float(sys.argv[2])
    avg_run, tau_t, dicts, cache, chemo_rates = run_sim(occupancies, dicts, chemo_rates, cum_chemo_rates, T, Dn=Dn, N=N, Dc=Dc, Nc=Nc, gamma=gamma, kappa=kappa, mu=mu, lambd=lambd, save_rate=1e5)

    save_sim(avg_run, tau_t, cache, chemo_rates)

if __name__ == '__main__':
    main()