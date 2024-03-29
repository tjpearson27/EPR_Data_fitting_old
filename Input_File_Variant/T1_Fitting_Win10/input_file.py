#function to use to fit the T1 data (monoexp, stretched, biexp, spec_diff)

function = 'stretched'

#general parameters (a and b should be 1 and 0 if normalized properly, 
#T1 is best guess at lowest temp)
a = 1
b = 0
T1 = 31000000
mono_bnds = [(0.6, 100, -0.2), (1.2, 10000000000, 0.5)]

###########################################################################################
#only for spec_diffus (spectral diffusion parameter)
q = 10
#sd_bnds = [(a_lower, T1_lower, q_lower, b_lower), (upper bounds same order)]
sd_bnds = [(0.6, 100, 10, -0.2), (1.2, 40000000, 100000000, 0.01)]

###########################################################################################
#only for stretched (stretch factor)
c = 0.5
#str_bnds = [(b_lower, a1_lower, T1_lower, c_lower), (upper bounds same order)]
str_bnds = [(), ()]

###########################################################################################
#params for biexp
a2 = 0.5
T1long = 30000000
T1short = 1000
#biexp_bnds = [(b_lower, a1_lower, T1_lower, c_lower), (upper bounds same order)]
biexp_bnds = [(), ()]

###########################################################################################
#folder to spit data out into
folder_name = 'T1_fit_{0}_4'.format(function)

#misc instructions
temps = [5, 10, 15, 20, 30, 40, 50, 60]
show_plots = False
show_temp_dep = True            
err_append = True               #appends esds to the end of the values and truncates to show error