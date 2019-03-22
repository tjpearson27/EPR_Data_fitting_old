
import Tkinter, tkFileDialog
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt 

filepath = 'G:\My Drive\Lab\UIUC EPR Data\Day 2\T1\data_fit_stretchedInvRec2\\temp_dep_stretched.txt'
if filepath == '':
    root = Tkinter.Tk()                      
    root.withdraw() 
    filepath = tkFileDialog.askopenfilename()

def J8(y):                                          #function definitions
    return integrate.quad(lambda x: (x**8)*
                (np.e**(x)/((np.e**(x)-1)**2)), 0, y)

def direct_process(x, a):
    return a * x

def raman_process(x, b, theta):
    return b * (((x/theta)**9) * J8(theta/x)[0])

def local_modes(x, c, delta):
    return c * np.e**(delta/x)/((np.e**(delta/x)-1)**2)

def fitfunc(T, a, b, c, d, e):              #a = dir, b = ram, c = loc, d = delta, e = theta
    y = direct_process(T, a) + raman_process(T, b, e) + local_modes(T, c, d)
    return y

vfitfunc = np.vectorize(fitfunc)            #scipy.integrate and scipy.optimize.curve_fit don't play nicely if it isn't vectorized
vraman = np.vectorize(raman_process)        
#vreg_fitfunc = np.vectorize(reg_fitfunc)

#################################################################
# parameters
#################################################################
direct = 10.5


raman = 1000000


local = 10**4


delta = 300


theta = 160

#################################################################
#
#################################################################

params = np.array([direct, raman, local, delta, theta])
#bnds = [(d_lower, r_lower, l_lower, dta_lower, t_lower),(d_upper, r_upper, l_upper, dta_upper, t_upper)]
a = np.loadtxt(filepath, delimiter=',')
#print a
T = a[:,0]
y = np.log10(1/(a[:,1]/10**9))
print a[:,1]
print y
reg_y = 1/a[:,1]
print reg_y
#print T, y

#popt, pcov = scipy.optimize.curve_fit(vfitfunc, T , y, params, bounds = bnds, gtol = 1e-15, max_nfev = 10000)
#perr = np.sqrt(np.diag(pcov))
#direct, raman, local, delta, theta = popt 
#d_err, r_err, l_err, dta_err, th_err = perr
plt.plot(T, y, 'go', T, np.log10(vfitfunc(T , params[0], params[1], params[2], params[3], params[4])), 'r--', 
                    T, np.log10(direct_process(T , params[0])), 'b--', 
                    T, np.log10(vraman(T , params[1], params[4])), 'c--',
                    T, np.log10(local_modes(T , params[2], params[3])), 'm--'
                    )

plt.gca().set_ylim(bottom=0, top=np.amax(y)+1)
plt.show()