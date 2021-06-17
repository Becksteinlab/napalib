from sklearn.neighbors import KernelDensity
import numpy as np

from tqdm import tqdm

def difference_distribution(data1, data2, n_points=100):
    
    x_range = np.linspace(min([data1.min(), data2.min()])-3,
                          max([data1.max(), data2.max()])+3,
                          n_points)
    
    d1 = x_range[1] - x_range[0]
    d2 = d1
    
    kde1 = KernelDensity(bandwidth=0.5, kernel='gaussian')
    kde2 = KernelDensity(bandwidth=0.5, kernel='gaussian')
    
    kde1.fit(data1[:,None])
    kde2.fit(data2[:,None])
    
    log1 = kde1.score_samples(x_range[:,None])
    log2 = kde2.score_samples(x_range[:,None])
    
    prob1 = np.exp(log1).astype(np.float32) #KDE of data1
    prob2 = np.exp(log2).astype(np.float32) #KDE of data2
    
    # note that both KDEs are already aligned 
    
    diff_prob = np.zeros_like(prob1)
    
    # outer product of reversed kdes so that diag with offset
    # can capture the desired lines
    mat = d1 * d2 * np.outer(prob1[::-1], prob2[::-1])
    
    offset = -np.argmin(np.abs(x_range))
    
    for i in range(diff_prob.shape[0]):
        diff_prob[i] = mat.diagonal(offset=offset).sum()
        offset += 1
    
    return x_range, diff_prob, prob1, prob2
                
