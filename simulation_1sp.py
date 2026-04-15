import numpy as np
from helper_functions import *
import sys

def main():
    # parameters
    L = 1
    Nsites = 100
    h = L / Nsites

    Dn1 = float(sys.argv[2])
    Dn2 = 0
    Dc = 1
    gamma1 = 2
    gamma2 = 0
    kappa = 1
    lambd1 = float(sys.argv[1])
    lambd2 = 0
    Dc = 1
    kappa = 1

        # homogeneous solution (unsteady)
    rhon1 = 2
    rhon2 = 2
    rhoc = rhon1 * gamma1 / kappa
    mu = Nsites / rhoc
    mu = 10

    Deff = (Dn1 - lambd1 * (T1(rhoc, mu) - T0(rhoc, mu)))
    print(f'Deff: {Deff}')
    mu, rhoc

    T = float(sys.argv[3])

    output_dir = f"/scratch.local/gtucci/micro/lambda1_{lambd1}/Dn1_{Dn1}/T_{T}"


    Dn1 = Dn1 / h ** 2
    Dn2 = Dn2 / h ** 2
    Dc = Dc / h ** 2
    lambd1 = lambd1 / h ** 2
    lambd2 = lambd2 / h ** 2

    occupancies_n1, dict_n1 = initialize_fields(Nsites, 0 / rhon1, rho=rhon1)
    occupancies_n2, dict_n2 = initialize_fields(Nsites, 0, rho=0)
    occupancies_c, dict_c = initialize_fields(Nsites, 0 / rhoc, rho=rhoc)
    chemo_rates1, cum_chemo_rates1 = initialize_chemo_rates(occupancies_n1, occupancies_c, mu, True)
    chemo_rates2, cum_chemo_rates2 = initialize_chemo_rates(occupancies_n2, occupancies_c, mu, False)

    N1, N2, Nc = int(np.sum(occupancies_n1)), int(np.sum(occupancies_n2)), int(np.sum(occupancies_c))

    occupancies = np.vstack((occupancies_n1, occupancies_n2, occupancies_c))
    chemo_rates = np.vstack((chemo_rates1, chemo_rates2))
    cum_chemo_rates = np.vstack((cum_chemo_rates1, cum_chemo_rates2))
    dicts = [dict_n1, dict_n2, dict_c]

    avg_run, tau_t, dicts, cache, chemo_rates = run_sim(occupancies, dicts, chemo_rates, 
                                                        cum_chemo_rates, T=T, Dn1=Dn1, Dn2=Dn2, N1=N1, N2=N2, Dc=Dc, Nc=Nc, 
                                                        gamma1=gamma1, gamma2=gamma2, kappa=kappa, mu=mu, lambd1=lambd1,lambd2=lambd2,
                                                        save_rate=1e5, save=True, output_dir=output_dir)
    
    # output_dir = f"../lambda_{lambd1}/D_{Dn1}/T_{T}"
    # save_sim(output_dir, avg_run, tau_t, cache, chemo_rates)

    print("Job finished")

if __name__ == '__main__':
    main()

    