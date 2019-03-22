Run in Python 2.7

Required files:
sigfig.py   (required for esd appending)
input_file.py
fit_T1_main.py

NOTE: data must be in the form of the normalized and phased .dat files as spit out from the matlab script that Matt wrote

you should also include the temps in the filename so that it fits them in the right order (the fitting uses the refined parameters from the first refinement as the inputs for the next temp)

you should only need to edit the file input_file.py

1) Edit input_file.py and put in necessary parameters and choose fit function (do this in Notepad++, atom, or vsc or something)
2) run fit_T1_main.py (you should use anaconda or jupyter or something that has all the common data science modules installed)
3) navigate to the data folder with the .dat files you want to fit
4) script should run and plots will pop up (if you selected the option)
5) data will be saved to .txt files that you can open in excel or origin in the directory that the scripts are running from
6) parameters will be saved with esds appended. if esd formatting is weird (ie 'value +- esd' instead of value(esd) that's because the error was at least one order of magnitude larger than the value of the parameter to which it corresponds
