from lolabrotation1.ferroptosis import model
from opt2q.noise import NoiseModel

sample_size = 10
# size of my heterogeneous cell population

# ------- Data -------
script_dir = os.path.dirname(__file__)
file_path = os.path.join(script_dir, 'SuiX_et_al_2018_HCT116_Ferroptosis_Data.csv')

cell_viability = pd.read_csv(file_path)
experimental_conditions = cell_viability[['Dose (uM)']]
len_ec = experimental_conditions.shape[0]
data = DataSet(data=cell_viability[['Dose (uM)', 'Time', 'Cell Viability']], # make me a dataset with cell viability
               measured_variables={'Cell Viability': 'nominal'}) # alive or dead (nominal)
# data.measurement_error_df = cell_viability[['', 'stdev']]
# removing this line forces the model to a default value -- 20% variability, presents cell viability in variances

# Parameter values (formatted for noise model)
kc0, kc2, kf3, kc3, kf4, kr7 = (1.0e-05, 1.0e-02, 3.0e-08, 1.0e-02, 1.0e-06, 1.0e-02)
k_values = pd.DataFrame([['kc0', kc0, True],
                         ['kc2', kc2, True],   # co-vary with kc3
                         ['kf3', kf3, False],
                         ['kc3', kc3, True],
                         ['kf4', kf4, False],
                         ['kr7', kr7, False]],
                        columns=['param', 'value', 'apply_noise'])\
    .iloc[np.repeat(range(6), len_ec)]         # Repeat for each of the experimental treatments

ligand = pd.DataFrame(cell_viability['Dose (uM)'].values, columns=['value'])
ligand['param'] = 'L_0'
ligand['apply_noise'] = False

param_m = pd.concat([k_values, ligand], sort=False, ignore_index=True)
param_m['Dose (uM)'] = np.tile(cell_viability['Dose (uM)'].values, 7)

kc2_cv, kc3_cv, kc2_kc3_cor = (0.2, 0.2, 0.25)
kc2_var, kc3_var, kc2_kc3_cov = ((kc2 * kc2_cv) ** 2, (kc3 * kc3_cv) ** 2, kc2 * kc2_cv * kc3 * kc3_cv * kc2_kc3_cor)
param_cov = pd.DataFrame([['kc2', 'kc2', kc2_var],
                          ['kc3', 'kc3', kc3_var],
                          ['kc2', 'kc3', kc2_kc3_cov]],  # Covariance between 'kf3' and kf4'
                         columns=['param_i', 'param_j', 'value'])

# ------- Noise Model -------
NoiseModel.default_sample_size = sample_size
noise = NoiseModel(param_mean=param_m, param_covariance=param_cov)
parameters = noise.run()

# ------- Simulate dynamics -------
sim = Simulator(model=model, param_values=parameters, solver='cupsoda')
results = sim.run(np.linspace(0, 3600*5, 100))

# ------- Measurement model -------
fk = FractionalKilling(simulation_result=results,
                       dataset=data,
                       measured_values={'Cell Viability': ['PUFA_PL_liperox_obs', 'time']},
                       observables=['PUFA_PL_liperox_obs'],
                       experimental_conditions=cell_viability[['Dose (uM)']],
                       time_dependent=False)

fk.process.remove_step('log_scale')
fk.process.add_step(('ddx', Scale(columns='PUFA_PL_liperox_obs', scale_fn=derivative)), 0)
fk.process.add_step(('at_max_t', ScaleGroups(groupby='simulation', scale_fn=where_max, **{'var': 'cPARP_obs'})), 1)
fk.process.add_step(('log10', Scale(columns='PUFA_PL_liperox_obs', scale_fn='log10')), 2)
fk.process.add_step(('polynomial', Scale(columns=['PUFA_PL_liperox_obs', 'Time'], scale_fn=polynomial_features, **{'degree': 2})),
                    'standardize')  # add after the 'standardize' step
fk.process.get_step('classifier').n_jobs = 1  # set number of multiprocessing jobs
fk.process.get_step('classifier').classifier_kwargs = {'verbose': 5}
fk.setup()

viability_coef = np.array([[-6.0, -0.08267159, -0.01332687, -0.02024716, 0.01*0.12825571]])
viability_intercept = np.array([8.00282563])
measurement_model_params_ = {'classifier__coefficients__viability__coef_': viability_coef,
                             'classifier__coefficients__viability__intercept_': viability_intercept}
fk.process.set_params(**measurement_model_params_)


# -------- likelihood function -----------
@objective_function(noise_model=noise, simulator=sim, measurement_model=fk, return_results=False, evals=0)
def likelihood_fn(x):
    kc0_ = 10 ** x[0]                           # :  [(-7,  -3),    float  kc0
    kc2_ = 10 ** x[1]                           # :   (-5,   1),    float  kc2
    kf3_ = 10 ** x[2]                           # :   (-11, -6),    float  kf3
    kc3_ = 10 ** x[3]                           # :   (-5,   1),    float  kc3
    kf4_ = 10 ** x[4]                           # :   (-10, -4),    float  kf4
    kr7_ = 10 ** x[5]                           # :   (-8,   4),    float  kr7

    kc2_cv_ = x[6]                              # :   (0, 1),       float  kc2_cv
    kc3_cv_ = x[7]                              # :   (0, 1),       float  kc3_cv
    kc2_kc3_cor_ = x[8]                         # :   (-1, 1)       float  kc2_kc3_cor

    viability_coef = np.array([[x[9],          # :  (-100, 100),   float
                                x[10],         # :  (-100, 100),   float
                                x[11],         # :  (-100, 100),   float
                                x[12],         # :  (-100, 100),   float
                                x[13]]])       # :  (-100, 100),   float
    viability_intercept = np.array([x[14]])    # :  (-10, 10)]     float

    # Each process selects one of the 4 gpu
    process_id = current_process().ident % 4

    k_val = pd.DataFrame([['kc0', kc0_, True],
                          ['kc2', kc2_, True],  # co-vary with kc3
                          ['kf3', kf3_, False],
                          ['kc3', kc3_, True],
                          ['kf4', kf4_, False],
                          ['kr7', kr7_, False]],
                         columns=['param', 'value', 'apply_noise']) \
        .iloc[np.repeat(range(6), len_ec)]  # Repeat for each of the experimental treatments

    lig = pd.DataFrame(cell_viability['Dose (uM)'].values, columns=['value'])
    lig['param'] = 'L_0'
    lig['apply_noise'] = False

    param_mean = pd.concat([k_val, lig], sort=False, ignore_index=True)
    param_mean['Dose (uM)'] = np.tile(cell_viability['Dose (uM)'].values, 7)

    kc2_var_, kc3_var_, kc2_kc3_cov_ = ((kc2_ * kc2_cv_) ** 2, (kc3_ * kc3_cv_) ** 2, kc2_ * kc2_cv_ * kc3_ * kc3_cv_ * kc2_kc3_cor_)
    param_covariance = pd.DataFrame([['kc2', 'kc2', kc2_var_],
                                     ['kc3', 'kc3', kc3_var_],
                                     ['kc2', 'kc3', kc2_kc3_cov_]],  # Covariance between 'kf3' and kf4'
                                    columns=['param_i', 'param_j', 'value'])

    measurement_model_params = {'classifier__coefficients__viability__coef_': viability_coef,
                                'classifier__coefficients__viability__intercept_': viability_intercept}

    likelihood_fn.noise_model.update_values(param_mean=param_mean,
                                            param_covariance=param_covariance)

    simulator_parameters = likelihood_fn.noise_model.run()
    likelihood_fn.simulator.param_values = simulator_parameters

    likelihood_fn.simulator.sim.gpu = [process_id]
    sim_results = likelihood_fn.simulator.run(np.linspace(0, 5000, 100))

    likelihood_fn.measurement_model.update_simulation_result(sim_results)
    likelihood_fn.measurement_model.process.set_params(**measurement_model_params)

    likelihood_fn.evals += 1
    l = likelihood_fn.measurement_model.likelihood()

    print(likelihood_fn.evals)
    print(x)
    print(l)
    return l

