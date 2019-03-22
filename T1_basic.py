import numpy as np 
import matplotlib.pyplot as plt
import Tkinter, tkFileDialog
import scipy.optimize 

def fitfunc(t, a, b, om, d, T_osc):
    return a*(1-b*np.cos(om*t + d)*np.exp(-t/T_osc))


root = Tkinter.Tk()                             
root.withdraw() 
file_path = tkFileDialog.askopenfilename() 

file_ = np.loadtxt(file_path, delimiter=',')


a = 1
b = 1
om = 1
d = 1
T_osc = 1
params = np.array([a, b, om, d, T_osc])


t = file_[:,0]
data = file_[:,1]
fit = file_[:,2]

osc = data - fit 

popt, pcov = scipy.optimize.curve_fit(fitfunc, t, osc, params)
a_fit, b_fit, om_fit, d_fit, T_osc_fit = popt 
fit_params = np.array([a_fit, b_fit, om_fit, d_fit, T_osc_fit])

print fit_params
print fit_params[2]/(2*np.pi)

plt.plot(t, osc, t, fitfunc(t, a_fit, b_fit, om_fit, d_fit, T_osc_fit))
plt.show() 
#plt.plot(t, osc)
#plt.show()

