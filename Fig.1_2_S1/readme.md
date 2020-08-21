The files in this directory recreate Figs. 1, 2 and S1 in our preprint of the paper "Probing the Hidden Sensitivity of Intrinsically Disordered Proteins to their Chemical Environment" (https://www.biorxiv.org/content/10.1101/2020.08.17.252478v1).

The folder contains the following .py files:
- hidden_sensitivity.py contains functions that read in and process our experimental data. It is referenced by the scripts that create the figures.
- hidden_sensitivity_fig_1a.py recreates Fig. 1A and writes hidden_sensitivity_fig_1a.png to the working directory. 
- hidden_sensitivity_fig_1b.py recreates Fig. 1B and writes hidden_sensitivity_fig_1b.png to the working directory. 
- hidden_sensitivity_fig_2a.py recreates Fig. 2A and writes hidden_sensitivity_fig_2a.png to the working directory. 
- hidden_sensitivity_fig_2b.py recreates Fig. 2B and writes hidden_sensitivity_fig_2b.png to the working directory. 
- hidden_sensitivity_fig_s1a.py recreates Fig. S1A and writes hidden_sensitivity_fig_s1a.png to the working directory. 
- hidden_sensitivity_fig_s1b.py recreates Fig. S1B and writes hidden_sensitivity_fig_s1b.png to the working directory. 

The .csv files contain our experimental data.

The files 1a1.png, 1a2.png and 1a3.png are used for Fig. 1A.

To run the programs from the spyder IPython console:
- Put all of the files in the same directory.
- Set that directory as the spyder working directory.
- To run each script, type: %run ./[filename]
- For example, to create Fig. 1A, type: %run ./hidden_sensitivity_fig_1a.py

Contact: David Moses dmoses5@ucmerced.edu
