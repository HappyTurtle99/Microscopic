import os
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.special import iv
import numpy as np
import time
import datetime

class RandomDict:
    def __init__(self):
        self.d = {}
        self.keys = []
        self.index = {}
        self.val_keys = {}
        self.max_key = 0

    def copy(self):
        new = RandomDict()
        new.d = self.d.copy()
        new.keys = self.keys.copy()
        new.index = self.index.copy()
        new.val_keys = {v: lst.copy() for v, lst in self.val_keys.items()}
        new.max_key = self.max_key
        return new

    def insert(self, key, value):
        if key not in self.d:
            self.index[key] = len(self.keys)
            self.keys.append(key)

        old_val = self.d.get(key)

        self.d[key] = value

        if old_val is not None:
            self.val_keys[old_val].remove(key)

        self.val_keys.setdefault(value, []).append(key)

        self.max_key = max(self.max_key, key)

    #like insert but produces a new key as well
    def insert_value(self, value):
        key = self.max_key + 1
        self.insert(key, value)

    def delete(self, key):
        if key not in self.d:
            return

        value = self.d[key]
        self.val_keys[value].remove(key)

        if not self.val_keys[value]:
            del self.val_keys[value]

        idx = self.index[key]
        last_key = self.keys[-1]

        self.keys[idx] = last_key
        self.index[last_key] = idx

        self.keys.pop()
        del self.index[key]
        del self.d[key]

    def get_random(self):
        key = random.choice(self.keys)
        return key, self.d[key]

    def random_key_from_val(self, value):
        keys_at_site = self.val_keys[value]
        key = random.choice(keys_at_site)
        return key

# Variables named d or dict are instances of RandomDict.
# is the particle index and index its position on the lattice.

# iv(n, x) is nth order modified bessel function evaluated at x.

def report_progress(tau, T, next_progress, start_time):
    progress = tau / T
    elapsed = time.time() - start_time

    while progress >= next_progress and next_progress <= 1.0:
        print(
            f"{int(next_progress * 100)}% completed "
            f"(tau = {tau:.4f}/{T:.4f}) "
            f"| elapsed: {elapsed:.2f} s"
        )
        next_progress += 0.05

    return next_progress

def f(c, mu=1):
    return max(0, np.tanh(mu * c))
    # if mu * c > 100: #avoids overflow
    #     return 1
    # return np.exp(mu * c) / (1 + np.exp(mu * c)) 

def p(gamma, m):
    return np.exp(-2 * gamma) * iv(m, 2 * gamma)

def T0(gamma,mu=1, m=70):
    return sum([f(i * m, mu=mu) * p(gamma, i) for i in range(-m, m + 1)])

def T1(gamma, mu=1, m=70):
    return sum([i * f(i * m, mu=mu) * p(gamma, i) for i in range(-m, m + 1)])

def get_new_chemo_rates(rates, occupancies_n, occupancies_c, index, mu, attractive):
    rates = rates.reshape((rates.shape[0] // 2, 2))

    # sites affected by a change at site i are i-1, i, i+1:
    l, r = occupancies_c[(index - 1) % occupancies_c.shape[0]] - occupancies_c[index], occupancies_c[(index + 1) % occupancies_c.shape[0]] - occupancies_c[index]

    if not attractive:
        l = -l 
        r = -r
        #opposite differences means particles attracted to negative gradient in concentration.

    #jump rate 0 if no particle
    if occupancies_n[index] == 0:
        rates[index] = np.array([0, 0])
    else:
        rates[index] = np.array([f(l, mu=mu), f(r, mu=mu)])

    if occupancies_n[(index - 1) % occupancies_c.shape[0]] == 0:
        rates[(index - 1) % rates.shape[0], 1] = 0
    else:
        rates[(index - 1) % rates.shape[0], 1] = f(-l, mu=mu)

    if occupancies_n[(index + 1) % occupancies_c.shape[0]] == 0:
        rates[(index + 1) % rates.shape[0], 0] = 0
    else:
        rates[(index + 1) % rates.shape[0], 0] = f(-r, mu=mu)

   # now we want to update the cumulative ones in O(1) time
    rates = rates.ravel()
    cum_rates = np.cumsum(rates)

    return rates, cum_rates

def get_new_chemo_rates_diff(rates, occupancies_n, occupancies_c, indexi, indexf, mu, attractive):
    # by convention we need indexf = indexi +- 1
    rates, cum_rates = get_new_chemo_rates(rates, occupancies_n, occupancies_c, indexi, mu, attractive)
    rates, cum_rates = get_new_chemo_rates(rates, occupancies_n, occupancies_c, indexf, mu, attractive)

    return rates, cum_rates

def evap(occupancies, d, p_index):
    d_out = d
    index = d.d[p_index]
    change = np.zeros_like(occupancies)
    change[index] = -1
    d_out.delete(p_index)
    dcache = np.array([-1, 0])
    return occupancies + change, d_out, dcache

# for clarity p_index pertainins to dict_n
def create(occupancies_c, dict_c, dict_n, p_index):
    index = dict_n.d[p_index]
    dict_c_out = dict_c
    change = np.zeros_like(occupancies_c)
    change[index] = 1
    dcache = np.array([1, 0])
    dict_c_out.insert_value(index)
    return occupancies_c + change, dict_c_out, dcache

def diff(occupancies, d, p_index):
    r = np.random.rand()
    index = d.d[p_index]

    if occupancies[index] == 0:
        print(index, d)
        plt.plot(occupancies)

        raise Exception('Improbable diff')

    if r <= 1 / 2:
        target_i = (index + 1) % occupancies.shape[0]
        dcache = np.array([0, 1])
    else:
        target_i = (index - 1) % occupancies.shape[0]
        dcache = np.array([0, -1])

    # update occupancy array and dictionary
    occupancies[index], occupancies[target_i] = occupancies[index] - 1, occupancies[target_i] + 1
    d.delete(p_index)
    d.insert(p_index, target_i)

    return occupancies, d, dcache

def chemo(occupancies, d, cum_chemo_rates):
    dcache = np.array([0, 0])

    index, direction = get_random_chemo_index(d, cum_chemo_rates, )
    target_i = None

    if direction == 0:
        target_i = (index - 1) % occupancies.shape[0]
        dcache = np.array([0, -1j])
    elif direction == 1:
        target_i= (index + 1) % occupancies.shape[0]
        dcache = np.array([0, 1j])

    # update occupancies
    occupancies[index], occupancies[target_i] = occupancies[index] - 1, occupancies[target_i] + 1
    
    # decide which particle it is we just moved
    p_index = d.random_key_from_val(index)
    d.delete(p_index)
    d.insert(p_index, target_i)

    return occupancies, d, dcache, index

#solves the issue that you can have nonzero jump rate from empty site for some choices of rate function.
def get_random_chemo_index(d, cum_chemo_rates):

    r = np.random.rand() * cum_chemo_rates[-1]
    idx = np.searchsorted(cum_chemo_rates, r)
    indices = np.unravel_index(idx, (cum_chemo_rates.shape[0] // 2, 2))
    index = indices[0]
    direction = indices[1]

    if index not in d.val_keys:
        print(index, d.val_keys, cum_chemo_rates)
        raise Exception('chemotactic jump from empty site.')

    return index, direction

def timestep(occupancies, dicts, chemo_rates, cum_chemo_rates, **kwargs):
    Dn1 = kwargs['Dn1']
    Dn2 = kwargs['Dn2']
    Dc = kwargs['Dc']
    N1 = kwargs['N1']
    N2 = kwargs['N2']
    Nc = kwargs['Nc']
    gamma1 = kwargs['gamma1']
    gamma2 = kwargs['gamma2']
    kappa = kwargs['kappa']
    mu = kwargs['mu']
    lambd1 = kwargs['lambd1']
    lambd2 = kwargs['lambd2']

    r1 = np.random.rand()
    r2 = np.random.rand()

    cum_chemo_rates1, cum_chemo_rates2 = cum_chemo_rates[0], cum_chemo_rates[1]
    total_chemo_rate1 = cum_chemo_rates1[-1]
    total_chemo_rate2 = cum_chemo_rates2[-1]
    total_rate = Dn1 * N1 + Dn2 * N2 + Dc * Nc + gamma1 * N1 + gamma2 * N2 + kappa * Nc + lambd1 * total_chemo_rate1 + lambd2 * total_chemo_rate2
    dtau = (1 / total_rate) * np.log(1 / r1)
    dcache = None

    # should probably neaten up this code make it into several functions.
    if r2 <= Dn1 * N1 / total_rate: # n1 diff
        p_index, index = dicts[0].get_random()
        occupancies[0], dicts[0], dcache = diff(occupancies[0], dicts[0], p_index)
        indexf = (index + 1) % occupancies[2].shape[0] if dcache[1] == 1 else (index - 1) % occupancies[2].shape[0]
        chemo_rates[0], cum_chemo_rates[0] = get_new_chemo_rates_diff(chemo_rates[0], occupancies[0], occupancies[2], index, indexf, mu, True)
    elif r2 <= (Dn1 * N1 + Dn2 * N2) / total_rate: # n2 diff
        p_index, index = dicts[1].get_random()
        occupancies[1], dicts[1], dcache = diff(occupancies[1], dicts[1], p_index)
        indexf = (index + 1) % occupancies[2].shape[0] if dcache[1] == 1 else (index - 1) % occupancies[2].shape[0]
        chemo_rates[1], cum_chemo_rates[1] = get_new_chemo_rates_diff(chemo_rates[1], occupancies[1], occupancies[2], index, indexf, mu, False)
    elif r2 <= (Dn1 * N1 + Dn2 * N2 + Dc * Nc) / total_rate: # c diff
        p_index, index = dicts[2].get_random()
        occupancies[2], dicts[2], dcache = diff(occupancies[2], dicts[2], p_index)
        indexf = (index + 1) % occupancies[2].shape[0] if dcache[1] == 1 else (index - 1) % occupancies[2].shape[0]
        chemo_rates[0], cum_chemo_rates[0] = get_new_chemo_rates_diff(chemo_rates[0], occupancies[0], occupancies[2], index, indexf, mu, True)
        chemo_rates[1], cum_chemo_rates[1] = get_new_chemo_rates_diff(chemo_rates[1], occupancies[1], occupancies[2], index, indexf, mu, False)
    elif r2 <= (Dn1 * N1 + Dn2 * N2 + Dc * Nc + gamma1 * N1) / total_rate: # c creation by n1
        p_index, index = dicts[0].get_random()
        occupancies[2], dicts[2], dcache = create(occupancies[2], dicts[2], dicts[0], p_index)
        chemo_rates[0], cum_chemo_rates[0] = get_new_chemo_rates(chemo_rates[0], occupancies[0], occupancies[2], index, mu, True)
        chemo_rates[1], cum_chemo_rates[1] = get_new_chemo_rates(chemo_rates[1], occupancies[1], occupancies[2], index, mu, False)
    elif r2 <= (Dn1 * N1 + Dn2 * N2 + Dc * Nc + gamma1 * N1 + gamma2 * N2) / total_rate: # c creation by n2
        p_index, index = dicts[1].get_random()
        occupancies[2], dicts[2], dcache = create(occupancies[2], dicts[2], dicts[1], p_index)
        chemo_rates[0], cum_chemo_rates[0] = get_new_chemo_rates(chemo_rates[0], occupancies[0], occupancies[2], index, mu, True)
        chemo_rates[1], cum_chemo_rates[1] = get_new_chemo_rates(chemo_rates[1], occupancies[1], occupancies[2], index, mu, False)
    elif r2 <= (Dn1 * N1 + Dn2 * N2 + Dc * Nc + gamma1 * N1 + gamma2 * N2 + kappa * Nc) / total_rate: # c evaporation
        p_index, index = dicts[2].get_random()
        occupancies[2], dicts[2], dcache = evap(occupancies[2], dicts[2], p_index)
        chemo_rates[0], cum_chemo_rates[0] = get_new_chemo_rates(chemo_rates[0], occupancies[0], occupancies[2], index, mu, True)
        chemo_rates[1], cum_chemo_rates[1] = get_new_chemo_rates(chemo_rates[1], occupancies[1], occupancies[2], index, mu, False)
    elif r2 <= (Dn1 * N1 + Dn2 * N2 + Dc * Nc + gamma1 * N1 + gamma2 * N2 + kappa * Nc + lambd1 * total_chemo_rate1) / total_rate: # chemotaxis of n1
        try:
            occupancies[0], dicts[0], dcache, index = chemo(occupancies[0], dicts[0], cum_chemo_rates[0])
            indexf = (index + 1) % occupancies[2].shape[0] if dcache[1] == 1j else (index - 1) % occupancies[2].shape[0]
            chemo_rates[0], cum_chemo_rates[0] = get_new_chemo_rates_diff(chemo_rates[0], occupancies[0], occupancies[2], index, indexf, mu, True)            
        except Exception:
            raise Exception('here')[1]
    elif r2 <= (Dn1 * N1 + Dn2 * N2 + Dc * Nc + gamma1 * N1 + gamma2 * N2 + kappa * Nc + lambd1 * total_chemo_rate1 + lambd2 * total_chemo_rate2) / total_rate: # chemotaxis of n2
        try:
            occupancies[1], dicts[1], dcache, index = chemo(occupancies[1], dicts[1], cum_chemo_rates[1])
            indexf = (index + 1) % occupancies[2].shape[0] if dcache[1] == 1j else (index - 1) % occupancies[2].shape[0]
            chemo_rates[1], cum_chemo_rates[1] = get_new_chemo_rates_diff(chemo_rates[1], occupancies[1], occupancies[2], index, indexf, mu, False)            
        except Exception:
            raise Exception('here')
    elif r2 <= 1:
        raise Exception(f'baaaad rates !! --Santi')
    
    if Nc > 700:
        print(Dn1 * N1, Dn2 * N2, Dc * Nc, gamma1 * N1, gamma2 * N2, kappa * Nc, lambd1 * total_chemo_rate1, lambd2 * total_chemo_rate2)
        raise KeyboardInterrupt

    return occupancies, dicts, chemo_rates, cum_chemo_rates, dtau, dcache

def initialize_fields(N, deltarho, rho=1, hom=False, k=1):
    if hom:
        occupancies = rho * np.ones(N)
    else:
        k = k * (2 * np.pi) / N
        occupancies = np.array([int(rho * (1 + deltarho * np.exp(i * k * 1j).real)) for i in range(N)])

    d = RandomDict()
    p_index = 0
    for i, e in enumerate(occupancies):
        for j in range(e):
            d.insert(p_index, i)
            p_index += 1

    return occupancies, d

def initialize_chemo_rates(occupancies_n, occupancies_c, mu, attractive):
    chemo_rates = np.zeros((occupancies_c.shape[0], 2))
    chemo_rates = chemo_rates.ravel()
    cum_chemo_rates = np.zeros_like(chemo_rates)
    for i in range(occupancies_c.shape[0]):
        chemo_rates, cum_chemo_rates = get_new_chemo_rates(chemo_rates, occupancies_n, occupancies_c, i, mu, attractive)
    return chemo_rates, cum_chemo_rates

def run_parallel(occupancies, dicts, chemo_rates, cum_chemo_rates, N_runs, T=1, **kwargs):
    parallel_runs = []
    min_lth = float('inf')

    for i in range(N_runs):
        new_run, tau_t = run_sim(occupancies, dicts, chemo_rates, cum_chemo_rates, T, **kwargs)
        parallel_runs.append(new_run)
        min_lth = min(min_lth, len(new_run))

    parallel_runs_arr = None

    parallel_runs_arr = np.stack([run[:min_lth] for run in parallel_runs])
    avg_occupancies = np.average(parallel_runs_arr, axis=0)
    return avg_occupancies, tau_t[:min_lth]

def run_sim(occupancies, dicts, chemo_rates, cum_chemo_rates, **kwargs):
    kwargs = kwargs.copy()

    occupancies_out = occupancies.copy() #[occupancy.copy() for occupancy in occupancies]
    dicts_out = list(d.copy() if d is not None else None for d in dicts)
    chemo_rates_out = chemo_rates.copy()
    cum_chemo_rates_out = cum_chemo_rates.copy()

    if 'tau_t' not in kwargs:
        tau = 0
        tau_t = [tau]
    else:
        tau = kwargs['tau_t'][-1]
        tau_t = list(kwargs['tau_t'])

    if 'occupancies_t' not in kwargs:
        occupancies_t = []
        occupancies_t.append(occupancies_out.copy())
    else:
        occupancies_t = list(kwargs['occupancies_t'])

    if 'cache' not in kwargs:
        cache = np.array([0j,0j])
    else:
        cache = kwargs['cache']

    counter = 0

    start_time = time.time()
    next_progress = 0.05  # initialize once
    
    # run
    while tau <= kwargs['T']:
        try:
            occupancies_out, dicts_out, chemo_rates_out, cum_chemo_rates_out, dtau, dcache = timestep(occupancies_out, dicts_out, chemo_rates_out, cum_chemo_rates_out, **kwargs)

            tau += dtau
            next_progress = report_progress(tau.real, kwargs['T'], next_progress, start_time)
            cache += dcache
            kwargs['Nc'] += dcache[0]
        
            if counter % int(kwargs['save_rate']) == 0:
                occupancies_t.append(occupancies_out.copy())
                tau_t.append(tau)
                
            counter += 1

        except KeyboardInterrupt as e:
            occupancies_t.append(occupancies_out.copy())
            tau_t.append(tau)
            occupancies_t = np.array(occupancies_t)
            tau_t = np.array(tau_t)
            if kwargs['save']:
                save_sim(occupancies_t, tau_t, cache, chemo_rates_out, **kwargs)
            print(f'Interrupted, occupancies_t shape: {occupancies_t.shape, counter, e}, Nc: {kwargs['Nc']}')
            return occupancies_t, tau_t, dicts_out, cache, chemo_rates_out
    
    occupancies_t.append(occupancies_out.copy())
    tau_t.append(tau)
    
    print(f'Number of timesteps completed: {counter}, Nc: {kwargs['Nc']}')
    occupancies_t = np.array(occupancies_t)
    tau_t = np.array(tau_t)

    if kwargs['save']:
        save_sim(occupancies_t, tau_t, cache, chemo_rates_out, **kwargs)
    return occupancies_t, tau_t, dicts_out, cache, chemo_rates_out

def meso_avg(occupancies, w):
    output = np.zeros(occupancies.shape[0])
    for i, _ in enumerate(occupancies):
        indices = [(i - w + j) % occupancies.shape[0] for j in range(2 * w + 1)]
        average = sum([occupancies[index] for index in indices])
        average = average / len(indices)
        output[i] = average
    return output

def save_sim(avg_run, tau_t, cache, chemo_rates, **kwargs):

    output_dir = kwargs['output_dir']
    # Create directory with timestamp
    if os.path.isdir(output_dir):
        output_dir = output_dir + '_new'
        
    os.makedirs(output_dir, exist_ok=True)

    # Save files with clean names
    np.save(os.path.join(output_dir, 'occupancies_t.npy'), avg_run)
    np.save(os.path.join(output_dir, 'tau_t.npy'), tau_t)
    np.save(os.path.join(output_dir, 'cache.npy'), cache)
    np.save(os.path.join(output_dir, 'chemo_rates.npy'), chemo_rates)
    
    print(f"Saved to: {output_dir}")

def load_sim(dir):
    names_in_order = ['occupancies_t.npy', 'tau_t.npy', 'cache.npy', 'chemo_rates.npy']
    return (np.load(os.path.join(dir, file)) for file in names_in_order)

def plot_results(avg_run, tau_t, w, start, end, Nlines=5):
    # Optional: Ensure Matplotlib uses its built-in TeX parser
    plt.rcParams.update({
        "text.usetex": False, # Usually False is safer unless you have a full TeX distro installed
        "mathtext.fontset": "cm" # This makes it look like standard LaTeX math
    })

    # mesoscopic average (width 2 * w)
    run = avg_run[int(start * avg_run.shape[0]):int(end * avg_run.shape[0]), :, :]
    Nlines = min(run.shape[0], Nlines)
    colors = plt.cm.viridis(np.linspace(0, 1, Nlines))

    #n1 plot
    plt.figure(figsize=(10, 6)) # Create a fresh figure
    for i in range(Nlines):
        index = int(len(run) * (i) * (1 / (Nlines - 1)))
        if index == len(run):
            index -= 1
        data_to_plot = meso_avg(run[index, 0, :], w) # Ensure it's 1D
        # data_to_plot = avg_run[index, 0, :]
        plt.plot(
            data_to_plot,
            color=colors[i],
            label=fr"$\tau$ = {tau_t[index]}",
            alpha=0.7 # Makes overlapping lines easier to see
        )

    plt.axis('tight')
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.title('density of n1')
    plt.show()

    # n2 plot
    plt.figure(figsize=(10, 6)) # Create a fresh figure
    for i in range(Nlines):
        index = int(len(run) * (i) * (1 / (Nlines - 1)))
        if index == len(run):
            index -= 1
        data_to_plot = meso_avg(run[index, 1, :], w) # Ensure it's 1D
        # data_to_plot = avg_run[index, 0, :]
        plt.plot(
            data_to_plot,
            color=colors[i],
            label=fr"$\tau$ = {tau_t[index]}",
            alpha=0.7 # Makes overlapping lines easier to see
        )

    plt.axis('tight')
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.title('density of n2')
    plt.show()

    colors = plt.cm.plasma(np.linspace(0, 1, Nlines))

    plt.figure(figsize=(10, 6)) # Create a fresh figure

    # Optional: Ensure Matplotlib uses its built-in TeX parser
    plt.rcParams.update({
        "text.usetex": False, # Usually False is safer unless you have a full TeX distro installed
        "mathtext.fontset": "cm" # This makes it look like standard LaTeX math
    })

    for i in range(Nlines):
        index = int(len(run) * (i) * (1 / (Nlines - 1)))
        if index == len(run):
            index -= 1
        data_to_plot = meso_avg(run[index, 2, :], w) # Ensure it's 1D
        # data_to_plot = avg_run[index, 0, :]
        plt.plot(
            data_to_plot,
            color=colors[i],
            label=fr"$\tau$ = {tau_t[index]}",
            alpha=0.7 # Makes overlapping lines easier to see
        )

    plt.axis('tight')
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.title('chemical')
    plt.show()