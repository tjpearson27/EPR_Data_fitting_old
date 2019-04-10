import numpy as np

#function to use to fit the T2 data (str_Full_ESEEM, str_no_ESEEM)

function = 'str_Full_ESEEM'

# STRETCHED FULL ESEEM parameters (str_no_ESEEM uses only a, T2, d, e)
# a is the overall scaling factor
a = 1
a_lower = -np.inf
a_upper = np.inf

# b is the modulation amplitude (Mod_Amp)
b = 1
b_lower = 0
b_upper = np.inf

# om is the oscillation frequency (Osc_Freq). 2pi * Osc_Freq should give the Larmor Freq of the nucleus
#   with which the spin is interacting
om = 4.032
om_lower = -np.inf
om_upper = 10

# c is the modulation phase (Mod_ph)
c = 2.62e-02
c_lower = -np.inf
c_upper = np.inf

# T_osc is the ESEEM decay time
T_osc = 1
T_osc_lower = 0
T_osc_upper = 5

# T2
T2 = 2
T2_lower = -np.inf
T2_upper = np.inf

# d is the stretch factor (Str_fac)
d = 1
d_lower = -np.inf
d_upper = np.inf

# e is a general offset term (Offset)
e = 0
e_lower = -np.inf
e_upper = np.inf

full_ESEEM_bnds = [(a_lower, b_lower, om_lower, c_lower, T_osc_lower, T2_lower, d_lower, e_lower), 
                    (a_upper, b_upper, om_upper, c_upper, T_osc_upper, T2_upper, d_upper, e_upper)]

no_ESEEM_bnds = [(a_lower, T2_lower, d_lower, e_lower), 
                (a_upper, T2_upper, d_upper, e_upper)]
###########################################################################################
#folder to spit data out into
folder_name = 'T2_fit_{0}_2'.format(function)

#misc instructions
show_plots = True
show_temp_dep = True            #not implemented yet
err_append = True               #appends esds to the end of the values and truncates to show error