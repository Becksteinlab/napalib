from napalib.saltbridge.toptools import Residue, Pair
import numpy as np
import xarray as xr

import json
import os


def get_interdomain_pairs(da):
    """From a dataarray, get the list of pairs for all salt bridges"""
    pos = da.pos
    neg = da.neg

    all_pairs = [Pair(Residue.from_label(p), Residue.from_label(n)) for p in pos.values for n in neg.values]
    inter_domain = list(filter(lambda x: x.interdomain, all_pairs))
    return inter_domain


def inter_domain_pairs_filter_min_distance(da, cutoff=10):
    """Get a reduced form of the results from get_interdomain_pairs where only pairs that came within
    a cutoff are included
    """
    inter_domain = get_interdomain_pairs(da)
    return list(filter(lambda x: x.get_time_series(da).min() < cutoff, inter_domain))


def exclusive_pairs(pl1, pl2):
    return list(filter(lambda x: ((x in pl1) + (x in pl2)) == 1, pl1 + pl2))


def common_pairs(pl1, pl2):
    common = []
    for p in (pl1 + pl2):
        if ((p in pl1) and (p in pl2)) and (not (p in common)):
            common.append(p)
    return common


def pairs_frame_distance(pairs, da, f1, f2, cutoff=10):
    distances = np.array([da[[f1, f2], :, :].loc[:, str(p.positive), str(p.negative)] for p in pairs])
    distances[distances > cutoff] = cutoff

    return np.sqrt(((distances[:, 0] - distances[:, 1]) ** 2).sum())


def construct_distance_array(pairs, da, cutoff=10, stride=1):
    data = np.array([da.loc[::stride, str(p.positive), str(p.negative)] for p in pairs]).T
    data[data > cutoff] = 10

    n = data.shape[0]

    result = np.empty((n, n), dtype=np.float32)

    for i in range(n):
        stripe = np.sqrt(((data[i] - data[i:]) ** 2).sum(axis=1))
        result[i, i:] = stripe
        result[i:, i] = stripe

    return result


def read_cluster_file(filename):
    with open(filename, 'r') as F:
        return json.load(F)


def save_cluster_results(filename, clusters, stride):
    indices = [int(i) for i in clusters.cluster_centers_indices_]
    labels = [int(i) for i in clusters.labels_]

    data = {"labels": labels, "indices": indices}

    try:
        if not os.path.exists(filename):
            full = {stride: data}
            with open(filename, 'w') as F:
                json.dump(full, F)
                return True

        with open(filename, 'r') as F:
            full = json.load(F)

        full.update({stride: data})

        with open(filename, 'w') as F:
            json.dump(full, F)
            return True
    except Exception as e:
        raise e


def cluster_AP(input_file, min_distance=3, max_distance=10, stride=10, save_file=None):
    from sklearn.cluster import AffinityPropagation

    if save_file:
        data = read_cluster_file(save_file)
        if str(stride) in list(data.keys):
            return data[str(stride)]

    da = xr.load_dataarray(input_file)
    pairs = inter_domain_pairs_filter_min_distance(da, min_distance)
    distance = construct_distance_array(pairs, da, cutoff=max_distance, stride=stride)
    cluster = AffinityPropagation(affinity="precomputed")
    cluster.fit(-distance)

    if save_file:
        save_cluster_results(save_file, cluster, stride)

    return cluster
