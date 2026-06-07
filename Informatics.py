import pandas as pd
import numpy as np
def softmax(x):
    exp_x = np.exp(x)
    return exp_x / np.sum(exp_x)

def markov_process(linkers):
    AA = []
    for i in range(len(linkers)):
        aa = list(linkers.iloc[i]['Sequnence'])
        AA = list(set(AA + aa))
    AA.sort()
    s = len(AA)
    Transition_probability = {i:np.zeros(s) for i in AA}
    ind = {}
    for i in range(len(AA)):
        ind[AA[i]] = i
    for i in range(len(linkers)):
        aa = list(linkers.iloc[i]['Sequnence'])
        for j in range(len(aa)-1):
            Transition_probability[aa[j]][ind[aa[j+1]]] += 1
    Transition_probability = {i:softmax(Transition_probability[i]) for i in AA}
    return pd.DataFrame(Transition_probability),AA