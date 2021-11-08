#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''

###############################################################################
#
#  File:            AMIGOSIII.py
#  Authors:         Morgan Shine and Chengxin Zhang
#  Creation Date:   2021-10-28
#
#  Nucleic Rama Plot
#      Calculate the pseudo-torsion angles eta, theta, eta', and theta'  
#      and determine sugar pucker for selected RNA/DNA in PyMOL
#      
#      Output:
#          eta_theta_plot.png        - 2D plot of eta and theta
#          eta_theta_prime_plot.png  - 2D plot of eta' and theta'
#          nucleic_worm_plot.png         - 3D plot of sequence, eta, and theta
#          eta_theta.csv             - comma separated text for eta, theta, 
#                                      and sugar pucker
#          eta_theta_prime.csv       - comma separated text for eta', theta', 
#                                      and sugar pucker
#
#  Motif Searching
#      Generate a nucleic acid worm database from an input directory
#
#      Output (for each chain of each PDB file in input directory):
#          name_ch_worm.csv       - comma separated text for nucleic acid worm of 
#                                   chain ch of PDB name
#                                     
#
#      Perform a nucleic acid worm search using a probe worm and a nucleic acid worm database
#
#      Output:
#          name_worm_search.txt   - text file for all comparisons between 
#                                   probe name and files in nucleic acid worm database
#
###############################################################################

'''

#Import lines

from __future__ import division
from __future__ import generators
from __future__ import print_function

import sys
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymol import cmd
import pickle
import time
import os

if int(sys.version[0])<=2:
    from Tkinter import *
    from Tkinter import Tk
    from Tkinter import messagebox
    from Tkinter import filedialog
else:
    from tkinter import *
    from tkinter import Tk
    from tkinter import messagebox
    from tkinter import filedialog

###############################################################################

class ETPlot:

    def __init__(self, selection=None, name=None, symbols='', state=-1):
        if selection is not None:
            self.start(selection)
                
    def start(self, sel):
        self.lock = 1

        #Define directory to save all output files
        print("Please select a folder to save all output files.")
        application_window = Tk()
        output_location = filedialog.askdirectory(parent=application_window, initialdir=os.getcwd(), title="Please select a folder to save all output files.")
            
        #Eta vs theta plot and csv file 
        et_plot_answer = messagebox.askyesno("Question", "Would you like to plot eta vs theta?")
        
        if et_plot_answer == True:
            #Generate eta vs theta plot in matplotlib
            ydata = []
            zdata = []
            sugar_marker = []
            
            et_colors_space = {'et_colors': []}
            et_color_tuples = []
            
            for (model, index), (eta, theta, sugar) in ETPlot.get_etatheta(self, sel).items():
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "et_colors.append(color)", space=et_colors_space)
                ydata.append(eta)
                zdata.append(theta)
                sugar_marker.append(sugar)
            
            ydata = np.array(ydata)
            zdata = np.array(zdata)        
                            
            for i in et_colors_space['et_colors']:
                tmpcolor = cmd.get_color_tuple(i)
                et_color_tuples.append(tmpcolor)
                
            fig1 = plt.figure()
            ax1 = plt.gca()

            #Add shading for helical region
            ax1.axvspan(150, 190, alpha=0.1, color='gray')
            ax1.axhspan(190, 260, alpha=0.1, color='gray')

            #Add data points
            for i in range(len(ydata)):
                if sugar_marker[i] == "C3'-endo":
                    ax1.scatter(ydata[i], zdata[i], color=et_color_tuples[i], edgecolor='black', marker="o")
                elif sugar_marker[i] == "C2'-endo":
                    ax1.scatter(ydata[i], zdata[i], color=et_color_tuples[i], edgecolor='black', marker="^")
                else: 
                    ax1.scatter(ydata[i], zdata[i], color=et_color_tuples[i], edgecolor='black', marker="s")

            #Add other plot features
            ax1.set_xlabel("Eta")
            ax1.set_ylabel("Theta")
            ax1.set_xticks(np.arange(0, 370, 45))
            ax1.set_yticks(np.arange(0, 370, 45))
            ax1.grid(which='both')
            plt.show()
            
            #Save plot as png
            try:
                plt.savefig(output_location + "/eta_theta_plot.png")
            except OSError:
                print("File could not be saved because output folder is not writable.")
                print("Please select a new output folder.")
                application_window = Tk()
                output_location = filedialog.askdirectory(parent=application_window, initialdir=os.getcwd(), title="Please select a folder to save all output files.")
                plt.savefig(output_location + "/eta_theta_plot.png")
            
            #Save csv file with sequence, eta, theta, and sugar pucker
            model_list = []
            model_index_list = []
            eta_list = []
            theta_list = []
            sugar_list = []
            
            chain_list_space = {'chain_list': []}
            resi_list_space = {'resi_list': []}
            resn_list_space = {'resn_list': []}
            
            for (model, index), (eta, theta, sugar) in ETPlot.get_etatheta(self, sel).items():
                model_list.append(model)
                model_index_list.append((model, index))
                eta_list.append(eta)
                theta_list.append(theta)
                sugar_list.append(sugar)
            
            for (model, index) in model_index_list:
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "chain_list.append(chain)", space=chain_list_space)
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "resi_list.append(resi)", space=resi_list_space)
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "resn_list.append(resn)", space=resn_list_space)

            df = {'Model': model_list, 'Chain': chain_list_space["chain_list"], 'Resn': resn_list_space["resn_list"], 'Resi': resi_list_space["resi_list"], 'Eta': eta_list, 'Theta': theta_list, 'Sugar_Pucker': sugar_list}
            df = pd.DataFrame(data=df, columns=['Model', 'Chain', 'Resn', 'Resi', 'Eta', 'Theta', 'Sugar_Pucker'])
            df.to_csv(output_location + "/eta_theta.csv")
        
        #Eta' vs Theta' plot and csv file
        et_p_plot_answer = messagebox.askyesno("Question", "Would you like to plot eta' vs theta'?")
        
        if et_p_plot_answer == True:
            #Generate eta' vs theta' plot in matplotlib
            ydata_p = []
            zdata_p = []
            sugar_marker_p = []
            
            et_colors_space_p = {'et_colors_p': []}
            et_color_tuples_p = []
            
            for (model, index), (eta_prime, theta_prime, sugar) in ETPlot.get_etatheta_prime(self, sel).items():
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "et_colors_p.append(color)", space=et_colors_space_p)
                ydata_p.append(eta_prime)
                zdata_p.append(theta_prime)
                sugar_marker_p.append(sugar)
            
            ydata_p = np.array(ydata_p)
            zdata_p = np.array(zdata_p) 
                            
            for i in et_colors_space_p['et_colors_p']:
                tmpcolor = cmd.get_color_tuple(i)
                et_color_tuples_p.append(tmpcolor)
    
            fig2 = plt.figure()
            ax2 = plt.gca()

            #Add shading for helical region
            ax2.axvspan(150, 190, alpha=0.1, color='gray')
            ax2.axhspan(190, 260, alpha=0.1, color='gray')

            #Add data points
            for i in range(len(ydata_p)):
                if sugar_marker_p[i] == "C3'-endo":
                    ax2.scatter(ydata_p[i], zdata_p[i], color=et_color_tuples_p[i], edgecolor='black', marker="o")
                elif sugar_marker_p[i] == "C2'-endo":
                    ax2.scatter(ydata_p[i], zdata_p[i], color=et_color_tuples_p[i], edgecolor='black', marker="^")
                else:
                    ax2.scatter(ydata_p[i], zdata_p[i], color=et_color_tuples_p[i], edgecolor='black', marker="s")

            #Add other plot features 
            ax2.set_xlabel("Eta'")
            ax2.set_ylabel("Theta'")
            ax2.set_xticks(np.arange(0, 370, 45))
            ax2.set_yticks(np.arange(0, 370, 45))
            ax2.grid(which='both')
            plt.show()
            
            #Save plot as png
            try:
                plt.savefig(output_location + "/eta_theta_prime_plot.png")
            except OSError:
                print("File could not be saved because output folder is not writable.")
                print("Please select a new output folder.")
                application_window = Tk()
                output_location = filedialog.askdirectory(parent=application_window, initialdir=os.getcwd(), title="Please select a folder to save all output files.")
                plt.savefig(output_location + "/eta_theta_prime_plot.png")
                
            #Save csv file with sequence, eta', theta', and sugar pucker
            model_list_p = []
            model_index_list_p = []
            eta_prime_list = []
            theta_prime_list = []
            sugar_list_p = []
                        
            chain_list_space_p = {'chain_list': []}
            resi_list_space_p = {'resi_list': []}
            resn_list_space_p = {'resn_list': []}
            
            for (model, index), (eta_prime, theta_prime, sugar) in ETPlot.get_etatheta_prime(self, sel).items():
                model_list_p.append(model)
                model_index_list_p.append((model, index))
                eta_prime_list.append(eta_prime)
                theta_prime_list.append(theta_prime)
                sugar_list_p.append(sugar)
                
            for (model, index) in model_index_list_p:
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "chain_list.append(chain)", space=chain_list_space_p)
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "resi_list.append(resi)", space=resi_list_space_p)
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "resn_list.append(resn)", space=resn_list_space_p)
                
            df_p = {'Model': model_list_p, 'Chain': chain_list_space_p["chain_list"], 'Resn': resn_list_space_p["resn_list"], 'Resi': resi_list_space_p["resi_list"], "Eta'": eta_prime_list, "Theta'": theta_prime_list, 'Sugar_Pucker': sugar_list_p}
            df_p = pd.DataFrame(data=df_p, columns=['Model', 'Chain', 'Resn', 'Resi', "Eta'", "Theta'", 'Sugar_Pucker'])
            df_p.to_csv(output_location + "/eta_theta_prime.csv")
            
        #Generate nucleic acid worm plot
        worm_plot_answer = messagebox.askyesno("Question", "Would you like to generate a nucleic acid worm plot?")

        if worm_plot_answer == True:
            xdata_space = {'xdata': []}
            ydata = []
            zdata = []
            sugar_marker = []
            
            et_colors_space = {'et_colors': []}
            et_color_tuples = []
            
            for (model, index), (eta, theta, sugar) in ETPlot.get_etatheta(self, sel).items():
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "xdata.append(resi)", space=xdata_space)
                ydata.append(eta)
                zdata.append(theta)
                sugar_marker.append(sugar)
            
            xdata = [int(i) for i in xdata_space["xdata"]]
            xdata = np.array(xdata)
            ydata = np.array(ydata)
            zdata = np.array(zdata)        
                        
            for (model, index), (eta, theta, sugar) in ETPlot.get_etatheta(self, sel).items():
                cmd.iterate("sel and model " + str(model) + " and index " + str(index), "et_colors.append(color)", space=et_colors_space)
                
            for i in et_colors_space['et_colors']:
                tmpcolor = cmd.get_color_tuple(i)
                et_color_tuples.append(tmpcolor)
    
            fig = plt.figure()
            try: # older version of matplotlib does not have 3d projection
                ax = plt.axes(projection='3d')
       
                for i in range((len(xdata)-1)):
                    if xdata[i] < xdata[i+1]:
                        #If the eta or theta angle falls in the helical region, color the segment blue
                        if ((ydata[i] > 150 and ydata[i] < 190) and (zdata[i] > 190 and zdata[i] < 260)):
                            ax.plot3D(xdata[i:i+2], ydata[i:i+2], zdata[i:i+2], color='mediumblue')
                        #If the eta or theta angle falls outside the helical region, color the segment red
                        else:
                            ax.plot3D(xdata[i:i+2], ydata[i:i+2], zdata[i:i+2], color='red')
                            ax.plot3D(xdata[i-1:i+1], ydata[i-1:i+1], zdata[i-1:i+1], color='red')
                        ax.plot(xdata[i:i+2], ydata[i:i+2], zs=0, zdir='z', color='silver')
                        ax.plot(xdata[i:i+2], zdata[i:i+2], zs=360, zdir='y', color='silver')
                
                ax.set_xlabel("Sequence", fontsize=10)
                ax.set_ylabel("Eta", fontsize=10)
                ax.set_zlabel("Theta", fontsize=10)
                ax.set_yticks(np.arange(0, 370, 45))
                ax.set_zticks(np.arange(0, 370, 45))
                ax.tick_params(axis='both', which='major', labelsize=8)
                ax.grid(False)
                plt.show()

                #Save plot as png
                try:
                    plt.savefig(output_location + "/nucleic_worm_plot.png", dpi=300)
                except OSError:
                    print("File could not be saved because output folder is not writable.")
                    print("Please select a new output folder.")
                    application_window = Tk()
                    output_location = filedialog.askdirectory(parent=application_window, initialdir=os.getcwd(), title="Please select a folder to save all output files.")
                    plt.savefig(output_location + "/nucleic_worm_plot.png")
                    
            except Exception as error:
                print("WARNING! Old version of matplotlib does not support 3d projection.")
                print(error)
            
        self.lock = 0
        
    def get_etatheta(self, sel):
        #Define variables for eta, theta, sugar dictionary
        eta_theta_dict = {}
        model_index = []
        eta_theta_vals = []
        
        #Determine what objects are in the selection
        stored_chain_C4p_space = {'stored_chain_C4p': []}        
        cmd.iterate("sel and name C4'", 'stored_chain_C4p.append(chain)', space=stored_chain_C4p_space)
        
        #Check if atoms are already selected
        if (len(stored_chain_C4p_space["stored_chain_C4p"])) == 0:
            cmd.select('(all)')
            print("The Selector-Error has been corrected by selecting all atoms.")
        
        object_list = cmd.get_object_list("sele")
        
        #Save tmp_object.pdb file for each object in object_list
        tmp_object_list = []
        for obj in object_list:
            cmd.save(f"tmp_{obj}.pdb", "sele and " + str(obj), -1, "")
            tmp_object_list.append(f"tmp_{obj}.pdb")
        
        #Save current working directory
        directory = os.getcwd()

        #Check for NaTorsion
        NaTorsion = os.path.join(os.path.abspath(os.path.dirname(__file__)), "NaTorsion")
        if os.name=="nt":
            NaTorsion+=".exe"
        if os.path.exists(NaTorsion):
            cond = 1
        else:
            os.chdir(os.path.abspath(os.path.dirname(__file__)))
            os.system("g++ -O3 NaTorsion.cpp -o NaTorsion")
            NaTorsion = os.path.join(os.path.abspath(os.path.dirname(__file__)), "NaTorsion")
            if os.path.exists(NaTorsion):
                cond = 1
                os.chdir(directory)
            else:
                cond = 3
                print("Error: NaTorsion is not in the same directory as AMIGOSIII.py.")

        #Calculate eta, theta, and sugar pucker for each object in tmp_object_list
        for obj_idx in range(0, len(tmp_object_list)):
            #Run NaTorsion for filename
            if cond == 1:
                input = os.popen(NaTorsion + " " + str(tmp_object_list[obj_idx])).read()
            else:
                break 

            input = list(input.split("\n"))
            
            for line in input:
                if line.startswith("N c resi") or len(line) == 0:
                    continue
                if line[66:73] != "-360.00" and line[74:81] != "-360.00":
                    #Add model and index to model_index 
                    chain = line[2]
                    rnum = int(line[4:8].strip())
                    
                    tmpmodel_index = cmd.index("sel and " + str(object_list[obj_idx]) + " and chain " + str(chain) + " and resi " + str(rnum))[0]
                    model_index.append(tmpmodel_index) 
                    
                    #Find eta, theta, and sugar pucker and add them to eta_theta_vals
                    tmpeta = float(line[66:73].strip())
                    if tmpeta < 1:
                        tmpeta = tmpeta + 360
                    eta = np.around(tmpeta, decimals=1)
                    
                    tmptheta = float(line[74:81].strip())
                    if tmptheta < 1:
                        tmptheta = tmptheta + 360
                    theta = np.around(tmptheta, decimals=1)
                    
                    sugar_torsion = float(line[26:33].strip())
                    if sugar_torsion > 0:
                        sugar = "C3'-endo"
                    elif sugar_torsion < 0:
                        sugar = "C2'-endo"
                        
                    eta_theta_vals.append((eta, theta, sugar))
               
        #Add all model_index and eta_theta_vals pairs to eta_theta_dict        
        for i in range(len(model_index)):
            eta_theta_dict[model_index[i]] = eta_theta_vals[i]
                    
        return(eta_theta_dict)
    
    def get_etatheta_prime(self, sel):
        #Define variables for eta, theta, sugar dictionary
        eta_theta_p_dict = {}
        model_index = []
        eta_theta_p_vals = []
        
        #Determine what objects are in the selection
        stored_chain_C4p_space = {'stored_chain_C4p': []}        
        cmd.iterate("sel and name C4'", 'stored_chain_C4p.append(chain)', space=stored_chain_C4p_space)
        
        #Check if atoms are already selected
        if (len(stored_chain_C4p_space["stored_chain_C4p"])) == 0:
            cmd.select('(all)')
            print("The Selector-Error has been corrected by selecting all atoms.")
        
        object_list = cmd.get_object_list("sele")
        
        #Save tmp_object.pdb file for each object in object_list
        tmp_object_list = []
        for obj in object_list:
            cmd.save(f"tmp_{obj}.pdb", "sele and " + str(obj), -1, "")
            tmp_object_list.append(f"tmp_{obj}.pdb")
        
        #Save current working directory
        directory = os.getcwd()

        #Check for NaTorsion
        NaTorsion = os.path.join(os.path.abspath(os.path.dirname(__file__)), "NaTorsion")
        if os.path.exists(NaTorsion):
            cond = 1
        else:
            os.chdir(os.path.abspath(os.path.dirname(__file__)))
            os.system("g++ -O3 NaTorsion.cpp -o NaTorsion")
            NaTorsion = os.path.join(os.path.abspath(os.path.dirname(__file__)), "NaTorsion")
            if os.path.exists(NaTorsion):
                cond = 1
                os.chdir(directory)
            else:
                cond = 3
                print("Error: NaTorsion is not in the same directory as AMIGOSIII.py.")

        #Calculate eta, theta, and sugar pucker for each object in tmp_object_list
        for obj_idx in range(0, len(tmp_object_list)):
            #Run NaTorsion for filename
            if cond == 1:
                input = os.popen(NaTorsion + " " + str(tmp_object_list[obj_idx])).read()
            else:
                break 

            input = list(input.split("\n"))
            
            for line in input:
                if line.startswith("N c resi") or len(line) == 0:
                    continue
                if line[66:73] != "-360.00" and line[74:81] != "-360.00":
                    #Add model and index to model_index 
                    chain = line[2]
                    rnum = int(line[4:8].strip())
                    
                    tmpmodel_index = cmd.index("sel and " + str(object_list[obj_idx]) + " and chain " + str(chain) + " and resi " + str(rnum))[0]
                    model_index.append(tmpmodel_index) 
                    
                    #Find eta', theta', and sugar pucker and add them to eta_theta_p_vals
                    tmpeta_p = float(line[82:89].strip())
                    if tmpeta_p < 1:
                        tmpeta_p = tmpeta_p + 360
                    eta_p = np.around(tmpeta_p, decimals=1)
                    
                    tmptheta_p = float(line[90:97].strip())
                    if tmptheta_p < 1:
                        tmptheta_p = tmptheta_p + 360
                    theta_p = np.around(tmptheta_p, decimals=1)
                    
                    sugar_torsion = float(line[26:33].strip())
                    if sugar_torsion > 0:
                        sugar = "C3'-endo"
                    elif sugar_torsion < 0:
                        sugar = "C2'-endo"
                        
                    eta_theta_p_vals.append((eta_p, theta_p, sugar))
               
        #Add all model_index and eta_theta_vals pairs to eta_theta_dict        
        for i in range(len(model_index)):
            eta_theta_p_dict[model_index[i]] = eta_theta_p_vals[i]
                    
        return(eta_theta_p_dict)
        
class RNAworm:
    
    def __init__(self, state=-1):
        self.start()
        
    def start(self):
        #Generate a nucleic acid worm database if needed 
        database_answer = messagebox.askyesno("Question", "Would you like to generate a nucleic acid worm database?")
        
        if database_answer:
            #Select input directory for generate_database
            print("Please select the directory you would like to use to generate a nucleic acid worm database.")
            application_window = Tk()
            directory = filedialog.askdirectory(parent=application_window, initialdir=os.getcwd(), title="Please select the directory you would like to use to generate a nucleic acid worm database.")
        
            start1 = time.time()
            RNAworm.generate_database(self, directory)
            end1 = time.time()
            print(f"Time to generate nucleic acid worm database: {end1-start1} s")
            
        #Perform a nucleic acid worm search
        worm_search_answer = messagebox.askyesno("Question", "Would you like to perform a nucleic acid worm search?")
        
        if worm_search_answer:
            #Select nucleic acid worm database directory
            print("Please select the directory containing the nucleic acid worm database.")
            application_window = Tk()
            nucleic_worm_database = filedialog.askdirectory(parent=application_window, initialdir=os.getcwd(), title="Please select the directory containing the nucleic acid worm database.")
            
            probe_format = messagebox.askquestion("Question", "Is the nucleic acid probe worm a PyMOL object (select no) or a local file (select yes)?")
            if probe_format:
                print("Please select the nucleic acid probe worm.")
                application_window = Tk()
                probe = filedialog.askopenfilename(parent=application_window, title="Please select the nucleic acid probe worm.")
            else:
                #Write probe csv file for PyMOL selection 
                cmd.select("sele extend 6") #required to calculate eta and theta for first and last residues in original selection
                model_list = []
                index_list = []
                eta_list = []
                theta_list = []
                chain_number_list = []
                
                chain_list_space = {'chain_list': []}
                resi_list_space = {'resi_list': []}
                resn_list_space = {'resn_list': []}
                
                for (model, index), (eta, theta, sugar) in ETPlot.get_etatheta(self, sel="sele").items():
                    model_list.append(model)
                    index_list.append(index)
                    eta_list.append(np.around(eta, decimals=1))
                    theta_list.append(np.around(theta, decimals=1))
                    chain_number_list.append(' ')
                
                for index in index_list:
                    cmd.iterate(("sel and index " + str(index)), "chain_list.append(chain)", space=chain_list_space)
                    cmd.iterate(("sel and index " + str(index)), "resi_list.append(resi)", space=resi_list_space)
                    cmd.iterate(("sel and index " + str(index)), "resn_list.append(resn)", space=resn_list_space)
                    
                df = {'PDB': model_list, 'Chain_Number': chain_number_list, 'Chain': chain_list_space["chain_list"], 'NT_Number': resi_list_space["resi_list"], 'NT_ID': resn_list_space["resn_list"], 'Eta': eta_list, 'Theta': theta_list}
                df = pd.DataFrame(data=df, columns=['PDB', 'Chain_Number', 'Chain', 'NT_Number', 'NT_ID', 'Eta', 'Theta'])
                df.to_csv(nucleic_worm_database + "/tmp_probe.csv")
                
                #Define probe
                probe = "tmp_probe.csv"
            
            start2 = time.time()
            RNAworm.worm_search(self, probe, nucleic_worm_database)
            end2 = time.time()
            print(f"Time to perform nucleic acid worm search: {end2-start2} s")

    def generate_database(self, directory):
        os.chdir(directory)
        entries = os.listdir(directory)
        
        #Check for NaTorsion
        NaTorsion = os.path.join(os.path.abspath(os.path.dirname(__file__)), "NaTorsion")
        if os.path.exists(NaTorsion):
            cond = 1
        else:
            os.chdir(os.path.abspath(os.path.dirname(__file__)))
            os.system("g++ -O3 NaTorsion.cpp -o NaTorsion")
            NaTorsion = os.path.join(os.path.abspath(os.path.dirname(__file__)), "NaTorsion")
            if os.path.exists(NaTorsion):
                cond = 1
                os.chdir(directory)
            else:
                cond = 3
                print("Error: NaTorsion is not in the same directory as AMIGOSIII.py.")

        for entry in entries:
            if os.path.splitext(entry)[1] == '.pdb': 
                filename = entry
                rname = []
                chain = []
                rnum = []
                eta = []
                theta = []
                
                #Run NaTorsion for filename
                if cond == 1:
                    input = os.popen(NaTorsion + " " + str(filename)).read()
                else: 
                    break

                input = list(input.split("\n"))
  
                for line in input:
                    if line.startswith("N c resi") or len(line) == 0:
                        continue
                    if line[66:73] != "-360.00" and line[74:81] != "-360.00":
                        rname.append(line[0].capitalize())
                        chain.append(line[2])
                        rnum.append(int(line[4:8].strip()))
                        tmpeta = float(line[66:73].strip())
                        if tmpeta < 1:
                            tmpeta = tmpeta + 360
                        eta.append(np.around(tmpeta, decimals=1))
                        tmptheta = float(line[74:81].strip())
                        if tmptheta < 1:
                            tmptheta = tmptheta + 360
                        theta.append(np.around(tmptheta, decimals=1))

                rname = np.array(rname, dtype=object)
                chain = np.array(chain, dtype=object)
                rnum = np.array(rnum, dtype=object)
                eta = np.array(eta, dtype=object)
                theta = np.array(theta, dtype=object)

                unique_chain = np.unique(chain)
                
                unique_chain_num = []
                x=1
                for i in unique_chain:
                    unique_chain_num.append(x)
                    x=x+1
                unique_chain_num = np.array(unique_chain_num, dtype=object)
                
                for ch in unique_chain:
                    ch_indices = np.array(chain == ch)
                    
                    unique_rnum = np.unique(rnum[ch_indices])
                    len_unique_rnum = len(unique_rnum)
                    
                    tmprname = rname[ch_indices]
                    tmpchain = chain[ch_indices]
                    tmprnum = rnum[ch_indices]
                    tmpeta = eta[ch_indices]
                    tmptheta = theta[ch_indices]
                    
                    tmpPDB = np.repeat(os.path.splitext(filename)[0], len_unique_rnum)
                    ch_num = int(unique_chain_num[unique_chain == ch])
                    tmpchain_num = np.repeat(ch_num, len_unique_rnum)

                    df = {'PDB': tmpPDB, 'Chain_Number': tmpchain_num, 'Chain': tmpchain, 'NT_Number': tmprnum, 'NT_ID': tmprname, 'Eta': tmpeta, 'Theta': tmptheta}
                    df = pd.DataFrame(data=df, columns=['PDB', 'Chain_Number', 'Chain', 'NT_Number', 'NT_ID', 'Eta', 'Theta'])
                    df.to_csv(str(filename) + '_' + str(ch_num) + '_worm.csv')
        
    def worm_search(self, probe, directory):
        #Move to input directory
        os.chdir(directory)
        
        #Read probe worm file and store as pandas dataframe
        probe_name = str(probe)
        probe = pd.read_csv(probe) 
        
        #Record all database files
        entries = os.listdir(directory)
        
        entry_dict = {}
        
        for entry in entries:
            if (os.path.splitext(entry))[1] == '.csv':
                name = pd.read_csv(entry)
                entry_dict[entry] = np.array((name["PDB"], name["Chain_Number"], name["Chain"], name["NT_Number"], name["NT_ID"], name["Eta"], name["Theta"]))
        
        pickled_entries_out = open('pickled_entries', 'wb')
        pickle.dump(entry_dict, pickled_entries_out)
        pickled_entries_out.close()
        pickled_entries_in = open('pickled_entries', 'rb')
        new_entry_dict = pickle.load(pickled_entries_in)
        pickled_entries_in.close()
        
        #Define lists for final output file
        PDB_final = []
        chain_num_final = []
        chain_final = []
        start_final = []
        end_final = []
        sequence_final = []
        average_final = []
        change_et_final = []
        
        #Store probe eta and theta values in numpy arrays
        probe_eta_vals = np.array(probe["Eta"])
        probe_theta_vals = np.array(probe["Theta"])
            
        #Define M as length of probe worm
        M = len(probe_eta_vals)
        
        #Loop through each nucleic acid worm .csv file in the input directory
        for entry in new_entry_dict:
            database_file = new_entry_dict[entry]
            database_file_eta_vals = database_file[5]
            
            #Define L as length of database_file
            L = len(database_file_eta_vals)
            
            #Check if the database_file is shorter than the probe worm 
            if M > L:
                continue
            
            #Define N
            N = L-M+1
            
            #Define database_file variables
            database_file = new_entry_dict[entry]
            database_file_NT_Number = database_file[3]
            database_file_NT_ID = database_file[4]
            database_file_eta_vals = database_file[5]
            database_file_theta_vals = database_file[6]
            
            #Create matrices for probe eta and theta values
            probe_eta_matrix = np.tile(probe_eta_vals, N).reshape(N, M)
            probe_theta_matrix = np.tile(probe_theta_vals, N).reshape(N, M)
            
            #Create empty matrices for database_file eta and theta values
            database_worm_eta_matrix = np.zeros((N, M))
            database_worm_theta_matrix = np.zeros((N, M))
            start = np.zeros(N, dtype=int)
            end = np.zeros(N, dtype=int)
            sequence = np.zeros(N, dtype=object)
            
            #Fill in database_worm matrices with eta and theta values 
            for i in range(0, N):
                database_worm_eta_matrix[i] = database_file_eta_vals[i:i+M]
                database_worm_theta_matrix[i] = database_file_theta_vals[i:i+M]
                start[i] = database_file_NT_Number[i]
                end[i] = database_file_NT_Number[i+M-1]
                tmpsequence = database_file_NT_ID[i:i+M]
                sequence[i] = ''.join(tmpsequence)
                
            #Calculate delta(eta, theta) 
            diff_eta = abs(probe_eta_matrix - database_worm_eta_matrix)
            idxe = diff_eta > 180
            diff_eta[idxe] = 360 - diff_eta[idxe]
            
            diff_theta = abs(probe_theta_matrix - database_worm_theta_matrix)
            idxt = diff_theta > 360
            diff_theta[idxt] = 360 - diff_theta[idxt]
            
            sqr_diff_eta = diff_eta ** 2
            sqr_diff_theta = diff_theta ** 2
            
            change_et = np.around(np.sqrt(sqr_diff_eta + sqr_diff_theta), decimals=2)
            
            average = np.around(np.average(change_et, axis=1), decimals=2)
            
            #Collect information from current entry for final output
            PDB = database_file[0][0]
            PDB_array = np.repeat(PDB, N)
            
            chain_num = database_file[1][0]
            chain_num_array = np.repeat(chain_num, N)
            
            chain = database_file[2][0]
            chain_array = np.repeat(chain, N)
            
            #Append information from current entry to final lists
            PDB_final.extend(PDB_array.tolist())
            chain_num_final.extend(chain_num_array.tolist())
            chain_final.extend(chain_array.tolist())
            start_final.extend(start.tolist())
            end_final.extend(end.tolist())
            sequence_final.extend(sequence.tolist())
            average_final.extend(average.tolist())
            change_et_final.extend(change_et.tolist())

        #Create dictionary for delta(eta, theta) values   
        change_et_final_array = np.array(change_et_final)
        change_dict = {}
        for x in range(0, len(probe_eta_vals)):
            change_dict["Change_{0}".format(x+1)] = change_et_final_array[:, x]
        
        change_keys = list(change_dict.keys())
        
        #Create dataframe for final output
        df = {'PDB': PDB_final, 'Chain_Number': chain_num_final, 'Chain': chain_final, 'Start': start_final, 'End': end_final, 'Sequence': sequence_final, 'Average': average_final}

        for y in range(0, len(probe_eta_vals)):
            df[change_keys[y]] = change_dict[change_keys[y]]
            
        col_names = list(df.keys())
        
        df = pd.DataFrame(data=df, columns=col_names)
        df = df.drop_duplicates()
        df = df.sort_values(by=["Average"])
        df = df.reset_index(drop=True)
        
        #Generate output file    
        name_file = str(probe_name) + "_worm_search.txt"
        with open(name_file, 'w') as f: df.to_string(f, col_space=5)   
                
def Nucleic_Rama(sel='(all)', name=None, symbols='aa', filename=None, state=-1):

    dyno = ETPlot(sel, name, symbols, int(state))
    if filename is not None:
        dyno.canvas.postscript(file=filename)

# Extend these commands
cmd.extend('eta_theta_plot', Nucleic_Rama)
cmd.auto_arg[0]['eta_theta_plot'] = cmd.auto_arg[0]['zoom']

# Add to plugin menu

def __init_plugin__(self):
    self.menuBar.addcascademenu('Plugin', 'AMIGOS III', 'Plot Tools', label='AMIGOS III Tools')
    self.menuBar.addmenuitem('AMIGOS III', 'command', 'Launch Nucleic Rama Plot', label='Nucleic Rama Plot',
                             command=lambda: ETPlot('(enabled)'))
    self.menuBar.addmenuitem('AMIGOS III', 'command', 'Launch Motif Searching', label='Motif Searching',
                             command=lambda: RNAworm('enabled)'))

# vi:expandtab:smarttab
