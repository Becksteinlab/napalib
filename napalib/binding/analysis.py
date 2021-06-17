from napalib.system.universe import NapAUniverse

import numpy as np

def _prot2int(protomer):
    return 0 if protomer.upper() == 'A' else 1

def binding_free_energy(dataset, protomer, cutoff=5):
    prot = _prot2int(protomer)
    distances = dataset['distance'][:,prot,:].values
    N_points = distances.shape[0]
    
    bound = distances[distances <= cutoff]
    unbound = distances[distances > cutoff]
    
    try:
        return np.log(unbound.shape[0] / bound.shape[0])
        
    except ZeroDivisionError:
        print(f"Divison by zero, found no bound states with a cutoff of {cutoff}")
        return None
        
def unique_sodium_indices(dataset, protomer):
    prot = _prot2int(protomer)
    indices = ds['sodium_index'][:,prot,:].values
    return np.unique(indices)
