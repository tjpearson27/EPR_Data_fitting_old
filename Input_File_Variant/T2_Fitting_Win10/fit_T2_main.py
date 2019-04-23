# [x] reads parameters and instructions (ie what function to use) from T2_input_file.py
# [x] imports data from .dat files
# [x] performs data fitting as prescribed by T2_input_file.py
# [x] error analysis
# [x] saves data to file prescribed by input_file.py
# [x] saves parameters to file

# yet to do:
# temp dependence

import Tkinter, tkFileDialog
import glob 
import os 
import numpy as np 
import scipy.optimize
import matplotlib.pyplot as plt
import types,re,string
import T2_input_file
from sigfig import *

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

#define fit functions

def str_Full_ESEEM(t, a, b, om, c, T_osc, T2, d, e):
    return a*(1-b*np.cos(om*t + c)*np.exp(-t/T_osc))*np.exp(-(2*t/T2)**d) + e

def str_no_ESEEM(t, a, T2, d, e):
    return a*np.exp(-(2*t/T2)**d) + e

#import data

root = Tkinter.Tk()                                                         #opens all .dat files in directory of the file pointed to
root.withdraw() 
file_path = tkFileDialog.askopenfilename()      
file_dir = os.path.dirname(file_path)                                       #trims filename from paths
allfiles = glob.glob(file_dir + '/*.dat')                                   #aggregates all .dat files in path           
allfiles.sort()

#define fit func

fitfunc = T2_input_file.function

#import parameters from input file, import data from files and fit

directory = T2_input_file.folder_name
createFolder('./{0}/'.format(directory))

if fitfunc == 'str_Full_ESEEM':
    a = T2_input_file.a
    b = T2_input_file.b
    om = T2_input_file.om
    c = T2_input_file.c 
    T_osc = T2_input_file.T_osc 
    T2 = T2_input_file.T2
    d = T2_input_file.d
    e = T2_input_file.e 
    
    params = np.array([a, b, om, c, T_osc, T2, d, e])
    data = [('A', 'Mod_Amp', 'Osc_Freq', 'Mod_ph', 'T_osc', 'T2', 'Str_fac', 'Offset')]

    for file_ in allfiles:
        try:
            a = np.loadtxt(file_, skiprows = 1)                                             #takes normalized and phased files from Matt's matlab script output
            t = a[:,0]
            y = a[:,1]
            popt, pcov = scipy.optimize.curve_fit(str_Full_ESEEM, t, y, params, bounds = T2_input_file.full_ESEEM_bnds)     #least squares refinement of parameters against data
            a, b, om, c, T_osc, T2, d, e = popt                                             #sets next temperature guess based on prev. temp refinement
            perr = np.sqrt(np.diag(pcov))                                                   #converts covariance to esd
            a_err, b_err, om_err, c_err, T_osc_err, T2_err, d_err, e_err = perr 
            #errors.append((a_err, b_err, om_err, c_err, T_osc_err, T2_err, d_err, e_err))        

            if T2_input_file.err_append == True:
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

            raw_data = [t, y, str_Full_ESEEM(t, a, b, om, c, T_osc, T2, d, e)]                     #writes ASCII file that can be imported to origin for data plotting
            #temp_dep.append([temp, T2, T2_err])
            os.chdir(os.path.abspath(os.curdir) + '\{0}'.format(directory))
            np.savetxt('Data_and_Fit_StrFullESEEM_{0}.txt'.format(file_.split('\\')[-1].replace('.dat', '')), 
                    np.transpose(raw_data), delimiter = ',', fmt='%s')        
            os.chdir('..')
            if T2_input_file.show_plots == True:
                plt.plot(t, y, 'go', t, str_Full_ESEEM(t,  a, b, om, c, T_osc, T2, d, e), 'r--')                    #plots data
                plt.title(file_.split('\\')[-1], fontsize=20)
                
                plt.show()
        except: break


#t, a, T2, d, e
if fitfunc == 'str_no_ESEEM':
    a = T2_input_file.a
    T2 = T2_input_file.T2
    d = T2_input_file.d
    e = T2_input_file.e 
    
    params = np.array([a, T2, d, e])
    data = [('A', 'T2', 'Str_fac', 'Offset')]

    for file_ in allfiles:
        
        try:
            a = np.loadtxt(file_, skiprows = 1)                                             #takes normalized and phased files from Matt's matlab script output
            t = a[:,0]
            y = a[:,1]
            popt, pcov = scipy.optimize.curve_fit(str_no_ESEEM, t, y, params, bounds = T2_input_file.no_ESEEM_bnds)     #least squares refinement of parameters against data
            a, T2, d, e = popt                                             #sets next temperature guess based on prev. temp refinement
            perr = np.sqrt(np.diag(pcov))                                                   #converts covariance to esd
            a_err, T2_err, d_err, e_err = perr 
            #errors.append((a_err, b_err, om_err, c_err, T_osc_err, T2_err, d_err, e_err))        
            print T2_input_file.err_append
            if T2_input_file.err_append == True:
                a_str = round_sig_error(a, a_err, 1, paren = True)                              #appends the esd to the end of the value to which it refers
                b_str = round_sig_error(b, b_err, 1, paren = True)
                om_str = round_sig_error(om, om_err, 1, paren = True)
                c_str = round_sig_error(c, c_err, 1, paren = True)
                T_osc_str = round_sig_error(T_osc, T_osc_err, 1, paren = True)
                T2_str = round_sig_error(T2, T2_err, 1, paren = True)       
                d_str = round_sig_error(d, d_err, 1, paren = True)
                e_str = round_sig_error(e, e_err, 1, paren = True)
                data.append((a_str, b_str, om_str, c_str, T_osc_str, T2_str, d_str, e_str))
                print "good to here"
            else: 
                data.append((a, T2, d, e))

            raw_data = [t, y, str_no_ESEEM(t, a, T2, d, e)]                     #writes ASCII file that can be imported to origin for data plotting
            #temp_dep.append([temp, T2, T2_err])
            os.chdir(os.path.abspath(os.curdir) + '\{0}'.format(directory))
            np.savetxt('Data_and_Fit_StrNoESEEM_{0}.txt'.format(file_.split('\\')[-1].replace('.dat', '')), 
                    np.transpose(raw_data), delimiter = ',', fmt='%s')        
            os.chdir('..')
            
            if T2_input_file.show_plots == True:
                plt.plot(t, y, 'go', t, str_no_ESEEM(t,  a, T2, d, e), 'r--')                    #plots data
                plt.title(file_.split('\\')[-1], fontsize=20)
                
                plt.show()
        except: break



os.chdir(os.path.abspath(os.curdir) + '\{0}'.format(directory))
np.savetxt('parameters_T2_{0}.txt'.format(fitfunc), data, delimiter = ',',fmt='%s')