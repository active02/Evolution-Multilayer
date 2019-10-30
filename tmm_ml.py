import scipy as sp
from tmm import (coh_tmm, unpolarized_RT, ellips, position_resolved, find_in_structure_with_inf)
import random
from deap import creator, base, tools, algorithms
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def generate_tmm_input_parameter_set_from_core_nkt_list(angle_list, wavelength_list, nkt_list_core):
    nk_list = [1.0]
    t_list = [sp.inf]
    for i in range(0,len(nkt_list_core)):
        nk_list.append(nkt_list_core[i][0])
        t_list.append(nkt_list_core[i][1])
    nk_list.append(1.0)
    t_list.append(sp.inf)
    nkt_list = sp.array([sp.array(nk_list), sp.array(t_list)])
    return [angle_list, wavelength_list, nkt_list]

def complex_nk_t_list_from_n_k_t_list(n_list, k_list, t_list):
    nk_list = [complex(n_list[i], k_list[i]) for i in range(0,len(n_list))]
    nk_list = sp.insert(nk_list,0,1.0)
    nk_list = sp.append(nk_list,1.0)
    t_list = sp.insert(t_list,0,sp.inf)
    t_list = sp.append(t_list,sp.inf)
    return nk_list, t_list

#Calculate spectrum
def reflection_spectrum_from_input_parameter_set(input_parameter_set):
    angle_list = input_parameter_set[0]
    wavelength_list = input_parameter_set[1]
    nkt_list = input_parameter_set[2]
    result_array = []
    for angle in angle_list:
        for wavelength in wavelength_list:
            result_array.append([angle, wavelength, unpolarized_RT(nkt_list[0],nkt_list[1],angle,wavelength)['R']])
    return sp.array(result_array)

def reflection_spectrum_from_nkt_list(wavelength_list, n_list, k_list, t_list):
    nk_list, t_list = complex_nk_t_list_from_n_k_t_list(n_list, k_list, t_list)
    result_array = [[i,unpolarized_RT(nk_list,t_list,0,i)['R']] for i in wavelength_list]
    return sp.array(result_array)

#Monte_Carlo_Method
def difference_between_two_spectrum(spectrum1, spectrum2):
    error_list = sp.array(spectrum1) - sp.array(spectrum2)
    return sp.sqrt(sp.mean(error_list**2)),

def difference_between_tmm_and_target_spectrum(angle, nkt_list, target_spectrum):
    wavelength_list = target_spectrum[0]
    target_amplitude_list = target_spectrum[1]
    tmm_amplitude_list = reflection_spectrum_from_input_parameter_set(generate_tmm_input_parameter_set_from_core_nkt_list([angle], wavelength_list, nkt_list))
    return difference_between_two_spectrum(tmm_amplitude_list[:,2], target_amplitude_list)

#Evolutionary Algorithm
def optimize_tmm_nk_list_using_evolutionary_algorithm(number_of_layers, nkt_range, number_of_populations, number_of_generations, wavelength_division, target_spectrum):
    angle = 0.0

    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    def random_nkt(n_range,k_range,t_range):
        return [complex(n_range[0]+(n_range[1]-n_range[0])*random.random(),k_range[0]+(k_range[1]-k_range[0])*random.random()),t_range[0]+(t_range[1]-t_range[0])*random.random()]
    toolbox.register("attr_nkt", random_nkt, nkt_range[0], nkt_range[1], nkt_range[2])
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_nkt, n=number_of_layers)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    def evalUR(nkt_list):
        return difference_between_tmm_and_target_spectrum(angle, nkt_list, target_spectrum)
    toolbox.register("evaluate", evalUR)
    toolbox.register("mate", tools.cxOnePoint)
    toolbox.register("mutate", tools.mutShuffleIndexes, indpb=0.05)
    toolbox.register("select", tools.selBest, k=int(0.9*number_of_populations))

    population = toolbox.population(n=number_of_populations)
    NGEN=number_of_generations

    generations = []
    for gen in range(NGEN):
        # Select and clone the next generation individuals
        offspring = map(toolbox.clone, toolbox.select(population))

        # Apply crossover and mutation on the offspring
        offspring = algorithms.varAnd(offspring, toolbox, cxpb=1.0, mutpb=0.02)

        # Evaluate the individuals with an invalid fitness
        fitnesses = toolbox.map(toolbox.evaluate, offspring)
        for fit, ind in zip(fitnesses, offspring):
            ind.fitness.values = fit

        # The population is entirely replaced by the offspring
        population[:] = offspring

        top1 = tools.selBest(population, k=1)
        top1_s = []
        for dna in top1:
            print(evalUR(dna)[0])
            dna_s = []
            for i in dna:
                dna_s.append([sp.real(i[0]),sp.imag(i[0]),sp.real(i[1])])
            top1_s.append(dna_s)
        generations.append(top1_s[0])
    return sp.array(generations)

def analyze(filename, target_spectrum):
    wavelength_list = target_spectrum[0]

    gens = pd.read_csv(filename)
    N_gens = gens.shape[1]-2

    n_dataset = sp.array(gens.loc[gens['minor']==0,'0':str(N_gens-1)])
    k_dataset = sp.array(gens.loc[gens['minor']==1,'0':str(N_gens-1)])
    t_dataset = sp.array(gens.loc[gens['minor']==2,'0':str(N_gens-1)])

    N_gens = t_dataset.shape[1]
    number_of_layers = t_dataset.shape[0]
    maximum_thickness = int(t_dataset.sum(axis=0).max())

    n_map = sp.ones((maximum_thickness,N_gens))
    k_map = sp.zeros((maximum_thickness,N_gens))
    for generation_id in range(0,N_gens):
        depth = 0
        for layer_id in range(0,number_of_layers):
            n_map[depth:depth+int(t_dataset[layer_id,generation_id]),generation_id] = n_dataset[layer_id,generation_id]
            k_map[depth:depth+int(t_dataset[layer_id,generation_id]),generation_id] = k_dataset[layer_id,generation_id]
            depth += int(t_dataset[layer_id,generation_id])

    fig1 = plt.figure(num=1, figsize=[10,7])

    fig1a = fig1.add_subplot(421)
    fig1b = fig1.add_subplot(422)
    fig1c = fig1.add_subplot(423)
    fig1d = fig1.add_subplot(424)
    fig1e = fig1.add_subplot(425)
    fig1f = fig1.add_subplot(426)
    fig1g = fig1.add_subplot(414)

    #fig.1(a-d) show structures of gens
    fig1a.imshow(n_map, extent=[0,150,0,50])
    fig1b.imshow(k_map, extent=[0,150,0,50])

    n_plot2_data = n_map[:,N_gens-1:N_gens].transpose()
    for i in n_plot2_data:
        fig1c.plot(sp.linspace(0,n_map.shape[0],n_map.shape[0]), i)

    k_plot2_data = k_map[:,N_gens-1:N_gens].transpose()
    for i in k_plot2_data:
        fig1d.plot(sp.linspace(0,k_map.shape[0],k_map.shape[0]), i)

    #fig.1(e-f) show reflection spectra of gens and comparision with target spectrum
    fig1e_generations_for_spectra = [0,1,-1]
    for i in fig1e_generations_for_spectra:
        reflection_spectra = reflection_spectrum_from_nkt_list(wavelength_list, n_dataset[:,i],k_dataset[:,i],t_dataset[:,i])
        fig1e.plot(wavelength_list, reflection_spectra[:,1], label=str(i))
    fig1e.plot(target_spectrum[0], target_spectrum[1], "r--", label="target")
    fig1e.legend(bbox_to_anchor=[1,1])

    fig1f_generation_for_spectra = sp.arange(0,N_gens,1)
    fig1f_reflecation_spectra = [reflection_spectrum_from_nkt_list(wavelength_list, n_dataset[:,i],k_dataset[:,i],t_dataset[:,i])[:,1] for i in fig1f_generation_for_spectra]
    print(fig1f_generation_for_spectra)
    fig1f.imshow(fig1f_reflecation_spectra, extent=[0,150,0,50])

    fom_list = [difference_between_two_spectrum(fig1f_reflecation_spectra[i],target_spectrum[1]) for i in range(0,len(fig1f_reflecation_spectra))]
    fig1g.plot(fig1f_generation_for_spectra, fom_list, "ko")

    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    filename = "test.csv"
    number_of_layers = 10
    wavelength_division = 10
    target_spectrum = [sp.linspace(380,780,wavelength_division),1.0*sp.ones(wavelength_division)]
    for i in range(2): target_spectrum[1][i] = 0.0

    nkt_range = [[1.0,4.0],[0.0,0.1],[5,2000]]
    number_of_populations = 100
    number_of_generations = 10
    generations = optimize_tmm_nk_list_using_evolutionary_algorithm(number_of_layers, nkt_range,  number_of_populations, number_of_generations, wavelength_division, target_spectrum)
    pd.Panel(generations).to_frame().to_csv(filename)

    analyze(filename, target_spectrum)
