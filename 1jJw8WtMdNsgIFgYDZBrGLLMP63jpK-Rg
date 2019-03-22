import Tkinter, tkFileDialog
import numpy as np 
import scipy.optimize
import scipy.integrate as integrate

import matplotlib.pyplot as plt
import types,re,string

epat = re.compile(r'^([^e]+)e(.+)$')

def round_sig(x, n):
   '''round floating point x to n significant figures'''
   if type(n) is not types.IntType:
      raise TypeError, "n must be an integer"
   try:
      x = float(x)
   except:
      raise TypeError, "x must be a floating point object"
   form = "%0." + str(n-1) + "e"
   st = form % x
   num,expo = epat.findall(st)[0]
   expo = int(expo)
   fs = string.split(num,'.')
   if len(fs) < 2:
      fs = [fs[0],""]
   if expo == 0:
      return num
   elif expo > 0:
      if len(fs[1]) < expo:
         fs[1] += "0"*(expo-len(fs[1]))
      st = fs[0]+fs[1][0:expo]
      if len(fs[1][expo:]) > 0:
         st += '.'+fs[1][expo:]
      return st
   else:
      expo = -expo
      if fs[0][0] == '-':
         fs[0] = fs[0][1:]
         sign = "-"
      else:
         sign = ""
      return sign+"0."+"0"*(expo-1)+fs[0]+fs[1]
      
def round_sig_error(x, ex, n, paren=False):
   '''Find ex rounded to n sig-figs and make the floating point x
   match the number of decimals.  If [paren], the string is
   returned as quantity(error) format'''
   stex = round_sig(ex,n)
   if stex.find('.') < 0:
      extra_zeros = len(stex) - n
      sigfigs = len(str(int(x))) - extra_zeros
      stx = round_sig(x,sigfigs)
   else:
      num_after_dec = len(string.split(stex,'.')[1])
      stx = ("%%.%df" % num_after_dec) % (x)
   if paren:
      if stex.find('.') >= 0:
         stex = stex[stex.find('.')+1:]
      return "%s(%s)" % (stx,stex.replace('0', ''))
   return stx,stex

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
    
    y = np.log10(direct_process(T, a) + raman_process(T, b, e) + local_modes(T, c, d))
    return y

def reg_fitfunc(T, a, b, c, d, e):             
    y = direct_process(T, a) + raman_process(T, b, e) + local_modes(T, c, d)
    return y

file_path = 'G:\My Drive\Lab\UIUC EPR Data\Day 2\T1\data_fit_specDif_invRec\\temp_dep_specDif.txt'
if file_path == '':
    root = Tkinter.Tk()                      
    root.withdraw() 
    file_path = tkFileDialog.askopenfilename()

vfitfunc = np.vectorize(fitfunc)            #scipy.integrate and scipy.optimize.curve_fit don't play nicely if it isn't vectorized
vraman = np.vectorize(raman_process)        
vreg_fitfunc = np.vectorize(reg_fitfunc)

#################################################################
# initial guesses
#################################################################
direct = 7
d_lower = 0.0
d_upper = 12

raman = 70000
r_lower = 1000
r_upper = 1500000

local = 10**7
l_lower = 10**4
l_upper = 10**9

delta = 300
dta_lower = 250
dta_upper= 450

theta = 75
t_lower = 60
t_upper = 200
#################################################################
#
#################################################################

params = np.array([direct, raman, local, delta, theta])
bnds = [(d_lower, r_lower, l_lower, dta_lower, t_lower),(d_upper, r_upper, l_upper, dta_upper, t_upper)]
a = np.loadtxt(file_path, delimiter=',')
T = a[:,0]
y = np.log10(1/(a[:,1]/10**9))
reg_y = 1/a[:,1]

popt, pcov = scipy.optimize.curve_fit(vfitfunc, T , y, params, bounds = bnds, gtol = 1e-15, max_nfev = 10000)
perr = np.sqrt(np.diag(pcov))
direct, raman, local, delta, theta = popt 
d_err, r_err, l_err, dta_err, th_err = perr

print 'Direct = {} +- {}'.format(direct, d_err)
print 'Raman = {} +- {}'.format(raman, r_err)
print 'Local = {} +- {}'.format(local, l_err)
print 'Delta = {} +- {}'.format(delta, dta_err)
print 'Theta = {} +- {}'.format(theta, th_err)

plt.plot(T, y, 'go', T, vfitfunc(T , popt[0], popt[1], popt[2], popt[3], popt[4]), 'r--', 
                    T, np.log10(direct_process(T , popt[0])), 'b--', 
                    T, np.log10(vraman(T , popt[1], popt[4])), 'c--',
                    T, np.log10(local_modes(T , popt[2], popt[3])), 'm--'
                    )

plt.gca().set_ylim(bottom=0, top=6)
plt.show()

d_str = round_sig_error(direct, d_err, 1, paren = True)                              
r_str = round_sig_error(raman, r_err, 1, paren = True)       
l_str = round_sig_error(local, l_err, 1, paren = True)
dta_str = round_sig_error(delta, dta_err, 1, paren = True)
th_str = round_sig_error(theta, th_err, 1, paren = True)

data = ['Direct_coeff = {0}'.format(d_str), 'Raman_coeff = {0}'.format(r_str), 'Local_mode_coeff = {0}'.format(l_str), 
                    'Local_mode_energy = {0}'.format(dta_str), 'Debye_temp = {0}'.format(th_str)]

data_and_fit = [T, y, vreg_fitfunc(T, direct, raman, local, delta, theta), direct_process(T, direct), 
                    vraman(T, raman, theta), local_modes(T, local, delta)]

data_and_fit_arr = np.asarray(data_and_fit)
data_and_fit_arr = np.transpose(data_and_fit_arr) 

#print data



np.savetxt('Temp_dep_parameters.txt', data, delimiter = ',', fmt='%s')
np.savetxt('Data_and_fits_temp_dep_for_plotting.txt', data_and_fit_arr, delimiter = ',', fmt='%s')