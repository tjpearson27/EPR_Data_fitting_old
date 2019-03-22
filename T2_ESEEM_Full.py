#=======================================================================================================================
# Fits Hahn Echo Decay curves including ESEEM modulation. Because of the number of parameters, it might be smart to run 
# the data through the T2_no_ESEEM and ESEEM_only scripts first to get an idea of starting points for each of these 
# parameters. 
#=======================================================================================================================

import Tkinter, tkFileDialog
import glob 
import os 
import numpy as np 
import scipy.optimize
import matplotlib.pyplot as plt 
import types,re,string

epat = re.compile(r'^([^e]+)e(.+)$')

root = Tkinter.Tk()                             #opens all .dat files in directory of the file pointed to
root.withdraw() 
file_path = tkFileDialog.askopenfilename()      
file_dir = os.path.dirname(file_path)           #trims filename from paths
allfiles = glob.glob(file_dir + '/*.dat')       #aggregates all .dat files in path           
allfiles.sort()

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

def fitfunc(t, a, b, om, c, T_osc, T2, d, e):
    return a*(1-b*np.cos(om*t + c)*np.exp(-t/T_osc))*np.exp(-(2*t/T2)**d) + e

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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# a*(1-b*np.cos(om*t + c)*np.exp(-t/T_osc))*np.exp(-(2*t/T2)**d) + e

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

# temps you collected at
temps = [5, 10, 15, 20, 30, 40, 50, 60, 80, 90]

# whether you want plots to pop up as data is fit. you'll have
#    to close one for the next dataset to be fit
show_plots = True
show_temp_dep = False
err_append = True
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bnds = [(a_lower, b_lower, om_lower, c_lower, T_osc_lower, T2_lower, d_lower, e_lower), 
        (a_upper, b_upper, om_upper, c_upper, T_osc_upper, T2_upper, d_upper, e_upper)]

temp_counter = 0
params = np.array([a, b, om, c, T_osc, T2, d, e])
data = [('A', 'Mod_Amp', 'Osc_Freq', 'Mod_ph', 'T_osc', 'T2', 'Str_fac', 'Offset')]
errors = [('a_err', 'b_err', 'om_err', 'c_err', 'T_osc_err', 'T2_err', 'd_err', 'e_err')]
temp_dep = []

createFolder('./data_fit_ESEEM_03042019/')

for file_ in allfiles:
    try:
      a = np.loadtxt(file_, skiprows = 1)                                             #takes normalized and phased files from Matt's matlab script output
      t = a[:,0]
      y = a[:,1]
      temp = temps[temp_counter]
      temp_counter += 1
      popt, pcov = scipy.optimize.curve_fit(fitfunc, t, y, params, bounds = bnds)     #least squares refinement of parameters against data
      a, b, om, c, T_osc, T2, d, e = popt                                             #sets next temperature guess based on prev. temp refinement
      perr = np.sqrt(np.diag(pcov))                                                   #converts covariance to esd
      a_err, b_err, om_err, c_err, T_osc_err, T2_err, d_err, e_err = perr 
      errors.append((a_err, b_err, om_err, c_err, T_osc_err, T2_err, d_err, e_err))        

      if err_append == True:
         a_str = round_sig_error(a, a_err, 1, paren = True)                              #appends the esd to the end of the value to which it refers
         b_str = round_sig_error(b, b_err, 1, paren = True)
         om_str = round_sig_error(om, om_err, 1, paren = True)
         c_str = round_sig_error(c, c_err, 1, paren = True)
         T_osc_str = round_sig_error(T_osc, T_osc_err, 1, paren = True)
         T2_str = round_sig_error(T2, T2_err, 1, paren = True)       
         d_str = round_sig_error(d, d_err, 1, paren = True)
         e_str = round_sig_error(e, e_err, 1, paren = True)
         data.append((a_str, b_str, om_str, c_str, T_osc_str, T2_str, d_str, e_str))
      else: 
         data.append((a, b, om, c, T_osc, T2, d, e))

      raw_data = [t, y, fitfunc(t, a, b, om, c, T_osc, T2, d, e)]                     #writes ASCII file that can be imported to origin for data plotting
      temp_dep.append([temp, T2, T2_err])
      os.chdir(os.path.abspath(os.curdir) + '\data_fit_ESEEM_03042019')
      np.savetxt('Data_and_Fit_{0}_ESEEM_Full.txt'.format(file_.split('\\')[-1].replace('.dat', '')), 
                  np.transpose(raw_data), delimiter = ',', fmt='%s')
      os.chdir('..')
      if show_plots == True:
         plt.plot(t, y, 'go', t, fitfunc(t,  a, b, om, c, T_osc, T2, d, e), 'r--')                    #plots data
         plt.title(file_.split('\\')[-1], fontsize=20)
         #plt.axis([8*10^-7, 1, -0.1, 1.1])

         plt.show()
    except: break


temp_dep = np.asarray(temp_dep)

if show_temp_dep == True:
    plt.semilogy(temp_dep[:,0], temp_dep[:,1], 'ro')
    plt.show()

os.chdir(os.path.abspath(os.curdir) + '\data_fit_ESEEM_03042019')
np.savetxt('parameters_{0}.txt'.format('T2_ESEEM_full'), data, delimiter = ',',fmt='%s')                      
np.savetxt('errors.txt', errors, delimiter = ',', fmt='%s')
np.savetxt('temp_dep_T2_ESEEM_Full_+err.txt', temp_dep, delimiter = ',', fmt='%s')