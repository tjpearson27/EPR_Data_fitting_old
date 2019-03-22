#===========================================================================================
# Reads in output .txt file from T2_no_ESEEM, subtracts the exponential decay, leaving just
# the modulation. Then fits the modulation to the ESEEM function. Also gives the Larmor 
# frequency of the nucleus giving rise to the ESEEM as output. 
#===========================================================================================

import numpy as np 
import matplotlib.pyplot as plt
import Tkinter, tkFileDialog
import scipy.optimize 
import os
import glob
import types,re,string
from scipy.signal import savgol_filter

epat = re.compile(r'^([^e]+)e(.+)$')

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

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

def fitfunc(t, a, b, om, d, T_osc):
    return a*(1-b*np.cos(om*t + d)*np.exp(-t/T_osc))


root = Tkinter.Tk()                             
root.withdraw() 
file_path = tkFileDialog.askopenfilename() 
file_dir = os.path.dirname(file_path)                   #trims filename from paths
allfiles = glob.glob(file_dir + '/*normph.txt')         #aggregates all .txt files in path           
allfiles.sort()

a = 0
b = 10
om = 5
d = -1.5
T_osc = 1.0
params = np.array([a, b, om, d, T_osc])

#=================================================================================
#Uncomment the following to just run one file
#=================================================================================
file_ = np.loadtxt(file_path, delimiter = ',')

t = file_[:,0]
data = file_[:,1]
fit = file_[:,2]

osc = data - fit 
osc_smooth = savgol_filter(osc, 13, 5) # window size 50, polynomial order 3


popt, pcov = scipy.optimize.curve_fit(fitfunc, t, osc_smooth, params)
a_fit, b_fit, om_fit, d_fit, T_osc_fit = popt 
fit_params = np.array([a_fit, b_fit, om_fit, d_fit, T_osc_fit])

print fit_params
print fit_params[2]/(2*np.pi)

plt.plot(t, osc_smooth, 'b--', 
        #t, osc, 'r-', 
        t, fitfunc(t, a_fit, b_fit, om_fit, d_fit, T_osc_fit)
        )
plt.show()
#=================================================================================
#
#=================================================================================

#=================================================================================
#Comment everything following if just running one file
#=================================================================================
#createFolder('./data_fit_ESEEM_only/')
#parameters = [('A', 'Mod_Amp', 'Mod_Freq', 'Mod_Ph', 'Decay_Time')]
#err_append = False 

#for f in allfiles:
#    file_ = np.loadtxt(f, delimiter=',')
#    t = file_[:,0]
#    data = file_[:,1]
#    fit = file_[:,2]
#
#    osc = data - fit 
#
#    popt, pcov = scipy.optimize.curve_fit(fitfunc, t, osc, params)
#    a, b, om, d, T_osc = popt
#    a_fit, b_fit, om_fit, d_fit, T_osc_fit = popt 
#    perr = np.sqrt(np.diag(pcov))
#    a_fit_err, b_fit_err, om_fit_err, d_fit_err, T_osc_fit_err = perr
#
#
#    fit_params = np.array([a_fit, b_fit, om_fit, d_fit, T_osc_fit])
#    fit_err = np.array([a_fit_err, b_fit_err, om_fit_err, d_fit_err, T_osc_fit_err])
#
#    if err_append == True:
#        a_str = round_sig_error(a, a_fit_err, 1, paren = True)                              #appends the esd to the end of the value to which it refers
#        b_str = round_sig_error(b, b_fit_err, 1, paren = True)
#        om_str = round_sig_error(om, om_fit_err, 1, paren = True)
#        d_str = round_sig_error(d, d_fit_err, 1, paren = True)
#        T_osc_str = round_sig_error(T_osc, T_osc_fit_err, 1, paren = True)
#        parameters.append((a_str, b_str, om_str, d_str, T_osc_str))
#    else:
#        parameters.append((a, b, om, d, T_osc))
#
#    raw_data = [t, osc, fitfunc(t, a, b, om, d, T_osc)]
#
    #print fit_params, fit_err
#    print 'ESEEM corresponds to Larmor frequency of approximately {} MHz'.format(fit_params[2]/(2*np.pi))
#
#    os.chdir(os.path.abspath(os.curdir) + '\data_fit_ESEEM_only')
#    np.savetxt('Data_and_Fit_{0}_ESEEM_only.txt'.format(f.split('\\')[-1].replace('.dat', '')), 
#                np.transpose(raw_data), delimiter = ',', fmt='%s')
#    os.chdir('..')
#
#    plt.plot(t, osc, t, fitfunc(t, a_fit, b_fit, om_fit, d_fit, T_osc_fit))
#    plt.show() 
#
#os.chdir(os.path.abspath(os.curdir) + '\data_fit_ESEEM_only')
#np.savetxt('parameters_ESEEM_only.txt', parameters, delimiter = ',', fmt = '%s')