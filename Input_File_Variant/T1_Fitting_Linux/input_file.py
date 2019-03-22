#function to use to fit the T1 data (monoexp, stretched, biexp, spec_diff)

function = 'spec_diff'

#general parameters (a and b should be 1 and 0 if normalized properly, 
#T1 is best guess at lowest temp)
a = 1
b = 0
T1 = 310000000

#only for spec_diffus (spectral diffusion parameter)
q = 10000000

#only for stretched (stretch factor)
c = 0.5

#params for biexp
a2 = 0.5
T1long = 30000000
T1short = 10000

#folder to spit data out into
folder_name = 'T1_fit_{0}'.format(function)

#misc instructions
show_plots = True
show_temp_dep = True            #not implemented yet
err_append = True               #appends esds to the end of the values and truncates to show error