#####################################################################################
#
#   TJP 7/18
#   for fitting normalized and phased inversion recovery curves to a biexponential 
#                       function
#
#   1) run Bruker files through Matt's MatLab script first
#   2) organize the files in a directory and rename by temp (005K.dat, 010K.dat, etc.)
#   3) run this script to fit to a monoexponential function
#           *The best guess you provide below should be the best guess for the lowest
#                temperature data. The guess updates for each temperature after 
#                refining on the previous temperature. 
#   4) the fit parameters for each temperature will be output to a file called 
#           'parameters_inv_rec.txt'
#   5) the data and the fit for plotting the curves and fits in origin will be 
#           output to files called 'Data_and_Fit_****.txt'
#
######################################################################################

import Tkinter, tkFileDialog
import glob 
import os 
import numpy as np 
import scipy.optimize
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

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
        
def fitfunc(t, b, a1, T1long, a2, T1short):                        #defines the fit function  
    return b-(a1*np.exp(-t/T1long) + a2*np.exp(-t/T1short))
                                 
root = Tkinter.Tk()                             #opens all .dat files in directory of the file pointed to
root.withdraw() 
file_path = tkFileDialog.askopenfilename()      
file_dir = os.path.dirname(file_path)           #trims filename from paths
allfiles = glob.glob(file_dir + '/*.dat')       #aggregates all .dat files in path           
allfiles.sort()

#==============================================================================================
#USER EDITS HERE
#==============================================================================================
b = 1                   #normalization factor best guess for lowest temp (usually 1) and bounds
b_lower = 0
b_upper = np.inf
a1 = 1                  
a1_lower = 0
a1_upper = np.inf
T1long = 10                   
T1long_lower = 0
T1long_upper = np.inf
a2 = 1000
a2_lower = 0
a2_upper = 100000
T1short = 1000000
T1short_lower = 0
T1short_upper = np.inf
show_plots = True     #if you want to see the plots of the fits
show_temp_dep = True 
temps = [5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90]
fit_temp_dep = False
err_append = False
#==============================================================================================
#Don't edit past here
#==============================================================================================

params = np.array([b, a1, T1long, a2, T1short])
bnds = [(b_lower, a1_lower, T1long_lower, a2_lower, T1short_lower), (b_upper, a1_upper, T1long_upper, a2_upper, T1short_upper)]

data = [('B', 'a1', 'T1_1', 'a2', 'T1_2')]
temp_counter = 0
temp_dep = []

createFolder('./data_fit_T1_biexp/')

for file_ in allfiles:
    
    a = np.loadtxt(file_, skiprows = 1)                                             #takes normalized and phased files from Matt's matlab script output
    t = a[:,0]
    y = a[:,1]
    temp = temps[temp_counter]
    temp_counter += 1

    popt, pcov = scipy.optimize.curve_fit(fitfunc, t, y, params, bounds = bnds)     #least squares refinement of parameters against data
    b, a1, T1long, a2, T1short = popt                                                                 #sets next temperature guess based on prev. temp refinement
    perr = np.sqrt(np.diag(pcov))                                                   #converts covariance to esd
    b_err, a1_err, T1_long_err, a2_err, T1_short_err = perr                 
    print '{} +- {}'.format(popt, perr)
    
    if show_plots == True:
        plt.plot(t, y, 'go', t, fitfunc(t, b, a1, T1long, a2, T1short), 'r--')                    #plots data
        plt.title(file_.split('\\')[-1], fontsize=20)
        #plt.axis([8*10^-7, 1, -0.1, 1.1])

        plt.show()
    
    if err_append == True:
        b_str = round_sig_error(b, b_err, 1, paren = True)                              #appends the esd to the end of the value to which it refers
        a1_str = round_sig_error(a1, a1_err, 1, paren = True)       
        T1_long_str = round_sig_error(T1long, T1_long_err, 1, paren = True)
        a2_str = round_sig_error(a2, a2_err, 1, paren = True)       
        T1_short_str = round_sig_error(T1short, T1_short_err, 1, paren = True)
        data.append((b_str, a1_str, T1_long_str, a2_str, T1_short_str))
    else: 
        data.append((b, a1, T1long, a2, T1short))

    raw_data = [t, y, fitfunc(t, b, a1, T1long, a2, T1short)]                                         #writes ASCII file that can be imported to origin for data plotting
    temp_dep.append([temp, T1long, T1short])
    os.chdir(os.path.abspath(os.curdir) + '\data_fit_T1_biexp')
    np.savetxt('Data_and_Fit_{0}.txt'.format(file_.split('\\')[-1].replace('.dat', '')), 
                np.transpose(raw_data), delimiter = ',', fmt='%s')
    os.chdir('..')
    


temp_dep = np.asarray(temp_dep)

if show_temp_dep == True:
    plt.semilogy(temp_dep[:,0], temp_dep[:,1], 'ro')
    plt.show()

os.chdir(os.path.abspath(os.curdir) + '\data_fit_T1_biexp')
np.savetxt('parameters_inv_rec.txt', data, delimiter = ',',fmt='%s')                      
np.savetxt('temp_dep.txt', temp_dep, delimiter = ',', fmt='%s')