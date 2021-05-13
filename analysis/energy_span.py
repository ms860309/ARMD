import numpy as np


def energy_span(energy_list):


    energy_profile = energy_list
    if energy_profile[0] != 0:
        print("Set initial state to reference state!")

    energy_ts = [ts for i,ts in enumerate(energy_profile) if i!=0 and i!=-1
                and energy_profile[i]>energy_profile[i-1] and energy_profile[i]>energy_profile[i+1] ]
    energy_int = [intermediate for intermediate in energy_profile if intermediate not in energy_ts]
    energy_rxn = energy_profile[-1]-energy_profile[0]
    energy_span_matrix = np.zeros((len(energy_int),len(energy_ts)))

    for j,ts in enumerate(energy_ts):
        for i,intermediate in enumerate(energy_int):
            if energy_profile.index(ts)<energy_profile.index(intermediate):
                energy_span_matrix[i,j]=ts-intermediate+energy_rxn
            else:
                energy_span_matrix[i,j]=ts-intermediate

    span = np.amax(energy_span_matrix)
    loc = np.where(energy_span_matrix==span)
    span_int = energy_int[loc[0][0]]
    span_ts = energy_ts[loc[1][0]]
    
    #TDTS and TDI are the two states that determine the energy span
    return f'Energy span = {span}, TDTS = {span_ts}, TDI = {span_int}, rxn energy = {energy_rxn}'


# path_1 = [0, 21.12, 12.81, 19.07, 5.32, 64.74, 21.21, 85.74, 29.05]
# path_2 = [1,-3,2,-4,5,-2]

# print(energy_span(path_1))
# print(energy_span(path_2))