# [x] reads parameters and instructions (ie what function to use) from input_file.py
# [x] imports data from .dat files
# [x] performs data fitting as prescribed by input_file.py
# [x] error analysis
# [x] saves data to file prescribed by input_file.py
# [x] saves parameters to file
# [x] temp dependence

# yet to do:

import Tkinter, tkFileDialog
import glob 
import os 
import numpy as np 
import scipy.optimize
import matplotlib.pyplot as plt
import types,re,string
import input_file
import sigfig 

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

#define fit functions

def monoexp(t, a1, T1, b):                        
    return -a1*(np.exp(-(t/T1)) - b - 1) 

def spec_diff(t, a1, T1, q, b):                        
    return -a1*(np.exp(-(t/T1 + np.sqrt(t/q))) - b - 1) 

def stretched(t, b, a1, T1, c):                         
    return b-a1*np.exp(-(t/T1)**c)

def biexp(t, b, a1, T1long, a2, T1short):                                   
    return b-(a1*np.exp(-t/T1long) + a2*np.exp(-t/T1short))

#import data

root = Tkinter.Tk()                                                         #opens all .dat files in directory of the file pointed to
root.withdraw() 
file_path = tkFileDialog.askopenfilename()      
file_dir = os.path.dirname(file_path)                                       #trims filename from paths
allfiles = glob.glob(file_dir + '/*.dat')                                   #aggregates all .dat files in path           
allfiles.sort()

#define fit func

fitfunc = input_file.function

#import parameters from input file, import data from files and fit

directory = input_file.folder_name
createFolder('./{0}/'.format(directory))
temps = input_file.temps
temp_counter = 0
temp_dep = []

if fitfunc == 'monoexp':
    a1 = input_file.a
    b = input_file.b
    T1 = input_file.T1
    params = np.array([a1, T1, b])
    data = [('a1', 'T1', 'b')]

    for file_ in allfiles:
        a = np.loadtxt(file_, skiprows = 1)                                  #takes normalized and phased files from Matt's matlab script output
        t = a[:,0]
        y = a[:,1]
        temp = temps[temp_counter]
        temp_counter += 1
        popt, pcov = scipy.optimize.curve_fit(monoexp, t, y, params, bounds = input_file.mono_bnds)         #least squares refinement of parameters against data
        a1, T1, b = popt          
        params = np.array([a1, T1, b])                                       #sets next temperature guess based on prev. temp refinement
        perr = np.sqrt(np.diag(pcov))                                        #converts covariance to esd
        a1_err, T1_err, b_err = perr
        
        if input_file.err_append == True:
            try:
                a1_str = sigfig.round_sig_error(a1, a1_err, 1, paren = True)       
            except:
                a1_str = '{0} +- {1}'.format(a1, a1_err)
            try:
                T1_str = sigfig.round_sig_error(T1, T1_err, 1, paren = True)
            except:
                T1_str = '{0} +- {1}'.format(T1, T1_err)
            try:
                b_str = sigfig.round_sig_error(b, b_err, 1, paren = True)
            except:    
                b_str = '{0} +- {1}'.format(b, b_err)
                
            data.append((a1_str, T1_str, b_str))                          #appends the esd to the end of the value to which it refers
        else:
            data.append((a1, T1, b))
        
        raw_data = [t, y, monoexp(t, a1, T1, b)]
        temp_dep.append([temp, T1, T1_err])
        os.chdir(os.path.abspath(os.curdir) + '\{0}'.format(directory))
        np.savetxt('Data_and_Fit_monoexp_{0}.txt'.format(file_.split('\\')[-1].replace('.dat', '')), 
                    np.transpose(raw_data), delimiter = ',', fmt='%s')
        
        
        print a1, T1, b
        
        plt.semilogx(t, y, 'go', t, monoexp(t, a1, T1, b), 'r--')            #plots data
        plt.title('{0} K {1}'.format(temp, fitfunc), fontsize=20)
        plt.savefig('InvRec_{0}.png'.format(file_.split('\\')[-1]), format = 'png')

        if input_file.show_plots == True:       
            plt.show()
        else:
            plt.close()
        
        os.chdir('..')
   

if fitfunc == 'spec_diff':
    a1 = input_file.a
    b = input_file.b
    T1 = input_file.T1
    q = input_file.q 
    params = np.array([a1, T1, q, b])
    data = [('a1', 'T1', 'q', 'b')]
    print "a1, T1, q, b"
    for file_ in allfiles:
        a = np.loadtxt(file_, skiprows = 1)                                            
        t = a[:,0]
        y = a[:,1]
        temp = temps[temp_counter]
        temp_counter += 1
        popt, pcov = scipy.optimize.curve_fit(spec_diff, t, y, params, bounds = input_file.sd_bnds)     
        a1, T1, q, b = popt
        params = np.array([a1, T1, q, b])                                                             
        perr = np.sqrt(np.diag(pcov))                                                  
        a1_err, T1_err, q_err, b_err = perr
        
        if input_file.err_append == True:
            try:
                a1_str = sigfig.round_sig_error(a1, a1_err, 1, paren = True)       
            except:
                a1_str = '{0} +- {1}'.format(a1, a1_err)
            try:
                T1_str = sigfig.round_sig_error(T1, T1_err, 1, paren = True)
            except:
                T1_str = '{0} +- {1}'.format(T1, T1_err)
            try:
                q_str = sigfig.round_sig_error(q, q_err, 1, paren = True)       
            except:
                q_str = '{0} +- {1}'.format(q, q_err)    
            try:
                b_str = sigfig.round_sig_error(b, b_err, 1, paren = True)
            except:    
                b_str = '{0} +- {1}'.format(b, b_err)
                
            data.append((a1_str, T1_str, q_str, b_str))                          #appends the esd to the end of the value to which it refers
        else:
            data.append((a1, T1, q, b))
        
        raw_data = [t, y, spec_diff(t, a1, T1, q, b)]
        temp_dep.append([temp, T1, T1_err])
        os.chdir(os.path.abspath(os.curdir) + '\{0}'.format(directory))
        np.savetxt('Data_and_Fit_spec_diff_{0}.txt'.format(file_.split('\\')[-1].replace('.dat', '')), 
                    np.transpose(raw_data), delimiter = ',', fmt='%s')        
        
        
        
        print a1, T1, q, b
        plt.semilogx(t, y, 'go', t, spec_diff(t, a1, T1, q, b), 'r--')                    
        plt.title('{0} K {1}'.format(temp, fitfunc), fontsize=20)
        plt.savefig('InvRec_{0}.png'.format(file_.split('\\')[-1]), format = 'png')

        if input_file.show_plots == True:
            plt.show()
        else:
            plt.close()
        os.chdir('..')
        


if fitfunc == 'stretched':
    a1 = input_file.a
    b = input_file.b
    T1 = input_file.T1
    c = input_file.c
    params = np.array([b, a1, T1, c])
    data = [('b', 'a1', 'T1', 'c')]
    
    for file_ in allfiles:
        a = np.loadtxt(file_, skiprows = 1)                                             
        t = a[:,0]
        y = a[:,1]
        temp = temps[temp_counter]
        temp_counter += 1
        popt, pcov = scipy.optimize.curve_fit(stretched, t, y, params, bounds = input_file.str_bnds)     
        b, a1, T1, c = popt
        params = np.array([b, a1, T1, c])                                                             
        perr = np.sqrt(np.diag(pcov))                                                   
        b_err, a1_err, T1_err, c_err = perr

        if input_file.err_append == True:
            try:
                a1_str = sigfig.round_sig_error(a1, a1_err, 1, paren = True)       
            except:
                a1_str = '{0} +- {1}'.format(a1, a1_err)
            try:
                T1_str = sigfig.round_sig_error(T1, T1_err, 1, paren = True)
            except:
                T1_str = '{0} +- {1}'.format(T1, T1_err)
            try:
                c_str = sigfig.round_sig_error(c, c_err, 1, paren = True)       
            except:
                c_str = '{0} +- {1}'.format(c, c_err)    
            try:
                b_str = sigfig.round_sig_error(b, b_err, 1, paren = True)
            except:    
                b_str = '{0} +- {1}'.format(b, b_err)
                
            data.append((b_str, a1_str, T1_str, c_str))                          #appends the esd to the end of the value to which it refers
        else:
            data.append((b, a1, T1, c))
        
        raw_data = [t, y, stretched(t, b, a1, T1, c)]
        temp_dep.append([temp, T1, T1_err])
        os.chdir(os.path.abspath(os.curdir) + '\{0}'.format(directory))
        np.savetxt('Data_and_Fit_stretched_{0}.txt'.format(file_.split('\\')[-1].replace('.dat', '')), 
                    np.transpose(raw_data), delimiter = ',', fmt='%s')
        
        print b, a1, T1, c
        
        plt.semilogx(t, y, 'go', t, stretched(t, b, a1, T1, c), 'r--')                    
        plt.title('{0} K {1}'.format(temp, fitfunc), fontsize=20)
        plt.savefig('InvRec_{0}.png'.format(file_.split('\\')[-1]), format = 'png')
        
        
        if input_file.show_plots == True:         
            plt.show()
        else: 
            plt.close()
        
        os.chdir('..')


if fitfunc == 'biexp':
    b = input_file.b
    a1 = input_file.a
    a2 = input_file.a2 
    T1long = input_file.T1long
    T1short = input_file.T1short
    params = np.array([b, a1, T1long, a2, T1short])
    data = [('b', 'a1', 'T1long', 'a2', 'T1short')]
    
    for file_ in allfiles:
        a = np.loadtxt(file_, skiprows = 1)                                            
        t = a[:,0]
        y = a[:,1]
        temp = temps[temp_counter]
        temp_counter += 1
        popt, pcov = scipy.optimize.curve_fit(biexp, t, y, params, bounds = input_file.biexp_bnds)     
        b, a1, T1long, a2, T1short = popt
        params = np.array([b, a1, T1long, a2, T1short])                                                            
        perr = np.sqrt(np.diag(pcov))                                                  
        b_err, a1_err, T1long_err, a2_err, T1short_err = perr
        
        if input_file.err_append == True:
            try:
                a1_str = sigfig.round_sig_error(a1, a1_err, 1, paren = True)       
            except:
                a1_str = '{0} +- {1}'.format(a1, a1_err)
            try:
                T1long_str = sigfig.round_sig_error(T1long, T1long_err, 1, paren = True)
            except:
                T1long_str = '{0} +- {1}'.format(T1long, T1long_err)
            try:
                T1short_str = sigfig.round_sig_error(T1short, T1short_err, 1, paren = True)
            except:
                T1short_str = '{0} +- {1}'.format(T1short, T1short_err)
            try:
                a2_str = sigfig.round_sig_error(a2, a2_err, 1, paren = True)       
            except:
                a2_str = '{0} +- {1}'.format(a2, a2_err)    
            try:
                b_str = sigfig.round_sig_error(b, b_err, 1, paren = True)
            except:    
                b_str = '{0} +- {1}'.format(b, b_err)
                
            data.append((b_str, a1_str, T1long_str, a2_str, T1short_str))                          #appends the esd to the end of the value to which it refers
        else:
            data.append((b, a1, T1long, a2, T1short))
        
        raw_data = [t, y, biexp(t, b, a1, T1long, a2, T1short)]
        temp_dep.append([temp, T1long, T1long_err, T1short, T1short_err])
        os.chdir(os.path.abspath(os.curdir) + '\{0}'.format(directory))
        np.savetxt('Data_and_Fit_biexp_{0}.txt'.format(file_.split('\\')[-1].replace('.dat', '')), 
                    np.transpose(raw_data), delimiter = ',', fmt='%s')        
       
        print b, a1, T1long, a2, T1short
        
        plt.semilogx(t, y, 'go', t, biexp(t, b, a1, T1long, a2, T1short), 'r--')                    
        plt.title('{0} K {1}'.format(temp, fitfunc), fontsize=20)
        plt.savefig('InvRec_{0}.png'.format(file_.split('\\')[-1]), format = 'png')
        
        if input_file.show_plots == True:          
            plt.show()
        else:
            plt.close()
        
        os.chdir('..')


temp_dep = np.asarray(temp_dep)
print temp_dep
os.chdir(os.path.abspath(os.curdir) + '\{0}'.format(directory))

if input_file.show_temp_dep == True:
    if input_file.function == 'biexp':
        plt.semilogy(temp_dep[:,0], temp_dep[:,1], 'ro', temp_dep[:,0], temp_dep[:,3], 'go')
    else: 
        plt.semilogy(temp_dep[:,0], temp_dep[:,1], 'ro')
    
    plt.savefig('Temp_dep_{0}.png'.format(fitfunc), format = 'png')
    plt.show()


np.savetxt('parameters_T1_{0}.txt'.format(fitfunc), data, delimiter = ',',fmt='%s')

if input_file.function == 'biexp':
    np.savetxt('temp_dep_{0}.txt'.format(fitfunc), temp_dep, delimiter = ',', fmt ='%s', 
                    header = 'Temp, T1, Error')
else:
    np.savetxt('temp_dep_{0}.txt'.format(fitfunc), temp_dep, delimiter = ',', fmt ='%s', 
                    header = 'Temp, T1_long, T1_long_Error, T1_short, T1_short_Error')