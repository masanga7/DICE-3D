# DICE coefficient per slice for nrrd images

import os
import nrrd      
import shutil
import statistics
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as cls
import ARR_functions
import Register_functions
from Register_functions import Path
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict
import colorama
from colorama import Fore, init
from skimage import measure, morphology
from skimage.feature import canny
from scipy import ndimage
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

global analysis, reference, master
master = tk.Tk()   
analysis = Path()
reference = Path()



def dates():
    
    
    class Path:
        def __init__(self):
            self.selected_directory = "None"
            self.folder_path = tk.StringVar()
            self.file_path = tk.StringVar()

        def browse_button(self):
            self.selected_directory = filedialog.askdirectory()
            self.folder_path.set(self.selected_directory)

        def browse_file(self):
            self.selected_file = filedialog.askopenfilename()
            self.file_path.set(self.selected_file)


    def change_phantom(event):
        global w4
        selection=w2.get()
        if selection == '3D':
            choices4 = ['all', 'right', 'center', 'left', 'right and center', 'right and left', 'left and center']
            w4 = ttk.Combobox(master, values=choices4)
            w4.grid(row=4, column=1)
            return w4.get()
      
        elif selection =='JASZCZAK':
            choices4 = ['1']
            w4 = ttk.Combobox(master, values=choices4)
            w4.grid(row=4, column=1)
            return w4.get()
        
    def saber(event):
        select=w1.get()
        if select == 'CT':
            misel='CT'
            return misel
        elif select == 'MR':
            misel='MR'
            return misel


    def run ():
        global analysis, reference, master, modality, activity, seg_desired, phantom, operator, date
        main_path_CT = reference.selected_directory
        main_path_PET = analysis.selected_directory
        modality = str(w1.get())
        phantom = str(w2.get())
        activity = str(w3.get())
        seg_desired =str(w4.get())
        operator =e4.get()
        date=e5.get()
        master.destroy()
        return activity, modality, seg_desired, phantom, operator



    
    master.title('DICE 3D')

    choices1 = ['CT', 'MR']
    w1 = ttk.Combobox(master, values=choices1)
    w1.grid(row=0, column=1)
    tk.Label(master, text="CT/MR").grid(row=0)
    w1.bind("<<ComboboxSelected>>", saber)

    choices2 = ['3D', 'JASZCZAK']
    w2 = ttk.Combobox(master, values=choices2)
    w2.grid(row=3, column=1)
    tk.Label(master, text="Phantom").grid(row=3)
    
    w2.bind("<<ComboboxSelected>>", change_phantom)
    tk.Label(master, text="Choose tubes/sphere").grid(row=4)    

    choices3 = ['yes', 'no']
    w3 = ttk.Combobox(master, values=choices3)
    w3.grid(row=5, column=1)
    tk.Label(master, text="Are all filled? (for one tube/sphere only)").grid(row=5)

    button1 = tk.Button(master, text="Select directory", command=reference.browse_button)
    lbl1 = tk.Label(master, textvariable=reference.folder_path)
    lbl1.grid(row=1, column=2)
    button1.grid(row=1, column=1)
    tk.Label(master, text="CT/MR directory").grid(row=1)

    button2 = tk.Button(master, text="Select directory", command=analysis.browse_button)
    lbl2 = tk.Label(master, textvariable=analysis.folder_path)
    lbl2.grid(row=2, column=2)
    button2.grid(row=2, column=1)
    tk.Label(master, text="PET directory").grid(row=2)

    
    e4 = tk.Entry(master)
    e4.grid(row=6, column=1)
    e4.insert(-1, 'John Doe')
    tk.Label(master, text="Operator").grid(row=6)

    e5 = tk.Entry(master)
    e5.grid(row=7, column=1)
    e5.insert(-1, '20/03/22' )
    tk.Label(master, text="date").grid(row=7)

    

    tk.Button(master, text='Run', command=run).grid(row=8, column=1, sticky=tk.W, pady=4)
    tk.mainloop()

    return  activity, modality, seg_desired, phantom, operator, date




#  Import data

activity, modality, seg_desired, phantom, operator, date = dates()
print ('las variables pasadas son:', activity, modality, seg_desired, phantom, operator, date)


main_path_CT = reference.selected_directory
main_path_PET = analysis.selected_directory



#  Create output folders and change directories

d = os.getcwd()                                                                     
p = r'Output_DICE'.format(d)                                                        
try:
    os.makedirs(p)                                                                   
except OSError:
    pass
path = os.path.join(d, p)                                                           

ptemp = r'additionalfiles_DICE'.format(d)                                           

try:
    os.makedirs(ptemp)                                                              
except OSError:
    pass
pathtemp = os.path.join(d, ptemp)                                                   

os.chdir(path)                                                                      


# Inputs

path_CT = os.path.join(main_path_CT, modality + '.nrrd')
path_PET = os.path.join(main_path_PET, 'pet.nrrd')                        


################# Importing images
#  Reads the .nrrd format files and saves the image and DICOM tags separately

image_CT = nrrd.read(path_CT)[0]                                                    
header_CT = nrrd.read(path_CT)[1]                                                   
name_image_CT = os.path.split(path_CT)[1]                                           
                                                                                    
image_PET = nrrd.read(path_PET)[0]                                                  
header_PET = nrrd.read(path_PET)[1]                                                 
name_image_PET = os.path.split(path_PET)[1]

   

################# PET threshold segmentation
# Applies a threshold to the PET image and generates the variables VOIs_PET, VOIs_REF, centroids and orders them by the centroid.

PET_voxel = float(header_PET['space directions'][2][2])
th_value_PET = 0.4                                                                 
print('Se ha aplicado al PET un umbral de',th_value_PET*100,'%', '¿Quieres otro? si o no')
info_umbral= str(input())

while (info_umbral != 'si' and info_umbral != 'no'):
    print ('por favor introduce "si" o "no"')
    info_umbral=str(input())

if info_umbral == 'si':
    print('¿Qué valor quieres (en %)?')
    umbral_usuario=float(input())
    th_value_PET=umbral_usuario/100

if info_umbral == 'no':
    th_value_PET = 0.4


if phantom == 'JASZCZAK':
    VOIs_PET, VOIs_REF, centroids = Register_functions.reference_seg(image_PET, th_value_PET,phantom)   

    
if phantom == '3D':
    VOIs_PET, VOIs_REF, centroids = Register_functions.reference_seg(image_PET, th_value_PET,phantom)   
    print (centroids)
    pairs = list(zip(VOIs_PET, VOIs_REF, centroids))                                                    
    pairs.sort(key=lambda x: x[2])                                                                      
    VOIs_PET, VOIs_REF, centroids = map(np.array, zip(*pairs))                                          



################### Create a PdfPages object
pdf = PdfPages('Nuevo - Register_DICE_index.pdf')                                                       
firstPage = plt.figure(figsize=(11.69, 8.27))
firstPage.clf()
title = 'DICE Co-register'
txt = 'DICE coefficient per slice fot the PET/' + modality + ' system.'
txt2 = 'Study made by ' + operator + ' on ' + date
firstPage.text(0.5, 0.9, title, transform=firstPage.transFigure, size=26, ha="center")
firstPage.text(0.05, 0.8, txt, transform=firstPage.transFigure, size=20, ha="left")
firstPage.text(0.05, 0.7, txt2, transform=firstPage.transFigure, size=20, ha="left")
pdf.savefig()



#################### CT threshold segmentation
# creates a vector "names" that contains the elements to analyze, for example: left, center and right or sphere Nº1

mod='water '+modality                                                                                      ###################################33
size_voxel = header_CT['space directions']

Threshold_CT = Register_functions.thresholdseg(image_CT, mod, phantom,size_voxel).astype('uint8')


os.chdir(pathtemp)                                                                                      
if phantom == 'JASZCZAK':
    names=['1','2','3','4','5','6']
else:
    names = ['left', 'center', 'right']                                                                 

nrrd.write('Th_' + modality + '.nrrd', Threshold_CT, header_CT, index_order='F')



############### Choose tube. input though terminal

if seg_desired == 'all':

    # Obtain segmentations for DICE
    for i in range(0, len(names)):                                                                                                                              
        os.chdir(pathtemp)                                                                                                                                      
        nrrd.write('Seg_ref_' + names[i] + '.nrrd', VOIs_REF[i], header_PET, index_order='F')                                                                   
        nrrd.write('VOI_PET_'+names[i]+'.nrrd', VOIs_PET[i], header_PET, index_order='F')
        nrrd.write('VOI_PET_' + names[i] + '_con_cabeceraPET.nrrd', VOIs_PET[i], header_PET, index_order='F')
        VOI_PET_cabecera = nrrd.read(os.path.join(pathtemp, 'VOI_PET_' + names[i] + '_con_cabeceraPET.nrrd'))[0]

        
        # Intersection TH+VOI
        ARR_functions.resampling(path_CT, os.path.join(pathtemp, 'Seg_ref_'+names[i]+'.nrrd'), pathtemp, 'VOI_REF_' + modality + '_' + names[i])                
        ref_CT = nrrd.read(os.path.join(pathtemp, 'VOI_REF_' + modality + '_'+names[i]+'_resampled.nrrd'))[0].astype('bool')                                    
        th_CT = np.copy(Threshold_CT)                                                                                                                           
        th_CT[~ref_CT] = 0                                                                                                                                      
        th_CT.astype('uint8')
        nrrd.write('VOI_' + modality + '_'+names[i]+'.nrrd', th_CT, header_CT, index_order='F')                                                                 

        # Intersection TH+VOI CT
        ARR_functions.resampling(path_PET, os.path.join(pathtemp, 'VOI_' + modality + '_'+names[i]+'.nrrd'), pathtemp, 'VOI_' + modality + '_'+names[i])

        # Intersection VOI_PET+VOI_CT
        VOI_CT = nrrd.read(os.path.join(pathtemp, 'VOI_' + modality + '_'+names[i]+'_resampled.nrrd'))[0].astype('bool')
        intersect = VOI_CT*VOIs_PET[i]
        os.chdir(path)
        nrrd.write('intersection_'+names[i]+'.nrrd', intersect, header_PET, index_order='F')

        # Number of slices calculation
        li_pet= 0
        li_ct=0

        for s in range(0, len(VOI_PET_cabecera[i][0, 0, :])):
            if np.mean(VOI_PET_cabecera[:, :, s]) == 0:
                li_pet = li_pet+1
            else:
                break
        
        for s in range(li_pet, len(VOI_PET_cabecera[i][0, 0, :])):
            if np.mean(VOI_PET_cabecera[:, :, s]) != 0 :      
                lf_pet = s
            else:
                break
        print('La imagen PET va desde el corte ', li_pet, 'hasta ', lf_pet)


        for r in range(0, len(VOI_CT[0, 0, :])):
            if np.mean(VOI_CT[:, :, r]) == 0:
                li_ct = li_ct + 1
            else:
                break

        for s in range(li_ct, len(VOI_CT[0, 0, :])):
            if np.mean(VOI_CT[:, :, s])!= 0 :
                lf_ct = s 
            else:
                break


        # Elección de los cortes inicial y final para calcular Dice
       
        if li_pet<li_ct:    
            li=li_pet
            inc_z=li_ct-li_pet
        else:               
            li=li_ct
            inc_z=li_pet-li_ct

        if lf_pet>lf_ct:
            lf=lf_pet
        else:
            lf=lf_ct



        # DICE calculation

        dice = []
        z = []

        for j in range(li, lf+1):       
            ct_slice = VOI_CT[:, :, j]
            pet_slice = VOI_PET_cabecera[i][:, :, j]
            intersect_slice = intersect[:, :, j]
            area_ct = np.sum(ct_slice == True)
            area_pet = np.sum(pet_slice == 1)   
            area_intersect = np.sum(intersect_slice == 1)   
            
            dicecoeff = (2 * area_intersect) / (area_ct + area_pet)
            dice.append(dicecoeff)
            z.append(j)    

        mean_dice = statistics.mean(dice)
        StandarDesviation_dice=statistics.stdev(dice)



        # Get ref file for calculate distance between CT and PET.
        Ref_CT=nrrd.read(os.path.join(pathtemp, 'VOI_' + modality + '_' + names[i] + '_resampled.nrrd'))[1]
        Ref_PET=nrrd.read(os.path.join(pathtemp, 'VOI_PET_' + names[i] + '_con_cabeceraPET.nrrd'))[1]
        

        # Load position left corner pixel
        PixelPosition_CT = (float(Ref_CT['space origin'][0]), float(Ref_CT['space origin'][1]), float(Ref_CT['space origin'][2]))
        PixelPosition_PET = (float(Ref_PET['space origin'][0]), float(Ref_PET['space origin'][1]), float(Ref_PET['space origin'][2]))
        ConstPixelSpacing_CT = (float(Ref_CT['space directions'][0][0]), float(Ref_CT['space directions'][1][1]), float(Ref_CT['space directions'][2][2]))
        ConstPixelSpacing_PET = (float(Ref_PET['space directions'][0][0]), float(Ref_PET['space directions'][1][1]), float(Ref_PET['space directions'][2][2]))


        # Calculo en mm del primer corte del CT y del PET y su diferencia.
       
        pos_PET=PixelPosition_PET[2]+ConstPixelSpacing_PET[2]*li_pet
        pos_CT=PixelPosition_CT[2]+ConstPixelSpacing_CT[2]*li_ct
        DIF_CORTES=abs(pos_PET-pos_CT)


        #print (dice)   
     
        # DICE's plot according to the slice
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)

        ax1.plot(z, dice, color=cls.to_rgba('C7', 0.7))                             # DICE
        ax1.axhline(y=mean_dice, color=cls.to_rgba('C9', 0.8))                      # DICE mean value
        ax1.axhline(y=0.93, color=cls.to_rgba('#FFA500', 0.8), linestyle='dotted')       
        ax1.fill_between(z, dice, color=cls.to_rgba('C7', 0.3))
        ax1.set_xlim([min(z), max(z)])                                              #Eje de abscisas del gráfico, desde zmin hasta zmax  
        ax1.set_ylim([0, 1])    
        ax1.set_xlabel('Slices', fontsize=15)
        ax1.set_ylabel('DICE', fontsize=15)
        ax1.set_title('Dice similarity index ' + names[i] + ' tube', fontsize=20)
        
        if round(mean_dice, 2) >= 0.93:   
            if DIF_CORTES>2 :
                color=('black','red', 'black')
                ax1.legend([f'DICE coefficient',f'Mean DICE= {round(mean_dice, 2)} \u00B1 {round(StandarDesviation_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm \nThe PET ' + modality + ' images are not well co-registered',
                        f'DICE threshold value: 0.93'],labelcolor=color,title=f'Co-registration properties')
            else:
                color=('black','green', 'black')
                ax1.legend([f'DICE coefficient', f'Mean DICE= {round(mean_dice, 2)} \u00B1 {round(StandarDesviation_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm \nThe PET ' + modality + ' images are well co-registered',
                        f'DICE threshold value: 0.93'], labelcolor=color, title=f'Co-registration properties')
        else:
            color=('black','red', 'black')
            ax1.legend([f'DICE coefficient']+ [f'Mean DICE= {round(mean_dice, 2)} \u00B1 {round(StandarDesviation_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm \nThe PET ' + modality + ' images are not properly co-registered',
                 f'DICE threshold value: 0.93'], labelcolor=color, title=f'Co-registration properties')        

        pdf.savefig(fig)
    pdf.close()
    print ('fin del programa')





if seg_desired != 'all':

    if activity == 'no':
        if phantom == 'JASZCZAK':
            if seg_desired == '1':
                names = names[0]
                VOIs_REF = VOIs_REF[0]          
                VOIs_PET = VOIs_PET [0]         
                centroids = centroids[0]
            elif seg_desired == '2':
                names = names[1]
                VOIs_REF = VOIs_REF[1]
                VOIs_PET = VOIs_PET [1]
                centroids = centroids[1]
            elif seg_desired == '3':
                names = names[2]
                VOIs_REF = VOIs_REF[2]
                VOIs_PET = VOIs_PET [2]
                centroids = centroids[2]
            elif seg_desired == '4':
                names = names[3]
                VOIs_REF = VOIs_REF[3]          
                VOIs_PET = VOIs_PET [3]         
                centroids = centroids[3]
            elif seg_desired == '5':
                names = names[4]
                VOIs_REF = VOIs_REF[4]
                VOIs_PET = VOIs_PET [4]
                centroids = centroids[4]
            elif seg_desired == '6':
                names = names[5]
                VOIs_REF = VOIs_REF[5]
                VOIs_PET = VOIs_PET [5]
                centroids = centroids[5]

        if phantom == '3D':
            if seg_desired == 'center' or 'right' or 'left':
                VOIs_REF = VOIs_REF[0]
                VOIs_PET = VOIs_PET [0]
                centroids = centroids[0]
            if seg_desired == 'center':
                names = names[1]
            elif seg_desired == 'right':
                names = names[2]
            elif seg_desired == 'left':
                names = names[0]

            

        nrrd.write('Seg_ref_' + names + '.nrrd', VOIs_REF, header_PET, index_order='F')
        nrrd.write('VOI_PET_' + names + '.nrrd', VOIs_PET, index_order='F')
        nrrd.write('VOI_PET_' + names + '_con_cabeceraPET.nrrd', VOIs_PET, header_PET, index_order='F')
        VOI_PET_cabecera = nrrd.read(os.path.join(pathtemp, 'VOI_PET_' + names + '_con_cabeceraPET.nrrd'))[0]
 
        
              
        # Intersection TH+VOI
        ARR_functions.resampling(path_CT, os.path.join(pathtemp, 'Seg_ref_' + names + '.nrrd'), pathtemp, 'VOI_REF_' + modality + '_' + names)

        # VOI CT resampling
        ref_CT = nrrd.read(os.path.join(pathtemp, 'VOI_REF_' + modality + '_' + names + '_resampled.nrrd'))[0].astype('bool')
        th_CT = np.copy(Threshold_CT)


        th_CT[~ref_CT] = 0
        th_CT.astype('uint8')
        nrrd.write('VOI_CT_' + names + '.nrrd', th_CT, header_CT, index_order='F')

        # Intersection TH+VOI CT
        ARR_functions.resampling(path_PET, os.path.join(pathtemp, 'VOI_' + modality + '_' + names + '.nrrd'), pathtemp,
                                 'VOI_' + modality + '_' + names)                                                                               # VOI CT resampling
        
        # Intersection VOI_PET+VOI_CT
        VOI_CT = nrrd.read(os.path.join(pathtemp, 'VOI_' + modality + '_' + names + '_resampled.nrrd'))[0].astype('bool')

        
        
        #intersect = VOI_CT * VOIs_PET
        intersect = VOI_CT * VOI_PET_cabecera
        nrrd.write('intersection_' + names + '.nrrd', intersect, header_PET, index_order='F')

        # Number of slices calculation
        li_pet= 0
        li_ct=0

       
        for s in range(0, len(VOI_PET_cabecera[0, 0, :])):
            if np.mean(VOI_PET_cabecera[:, :, s]) == 0:
                li_pet = li_pet+1
            else:
                break
        
        for s in range(li_pet, len(VOI_PET_cabecera[0, 0, :])):
            if np.mean(VOI_PET_cabecera[:, :, s]) != 0 :      
                lf_pet = s
            else:
                break
        print('La imagen PET va desde el corte ', li_pet, 'hasta ', lf_pet)
        
       
        

        for r in range(0, len(VOI_CT[0, 0, :])):
            if np.mean(VOI_CT[:, :, r]) == 0:
                li_ct = li_ct + 1
            else:
                break


        for s in range(li_ct, len(VOI_CT[0, 0, :])):
            if np.mean(VOI_CT[:, :, s])!= 0 :
                lf_ct = s 
            else:
                break

        
        # Elección de los cortes inicial y final para calcular Dice
       
        if li_pet<li_ct:    
            li=li_pet
            inc_z=li_ct-li_pet
        else:               
            li=li_ct
            inc_z=li_pet-li_ct

        if lf_pet>lf_ct:
            lf=lf_pet
        else:
            lf=lf_ct


        
        # DICE calculation
        dice = []
        z = []

        for j in range(li, lf+1):       
            ct_slice = VOI_CT[:, :, j]
            pet_slice = VOI_PET_cabecera[:, :, j]
            intersect_slice = intersect[:, :, j]
            area_ct = np.sum(ct_slice == True)
            area_pet = np.sum(pet_slice == 1)   
            area_intersect = np.sum(intersect_slice == 1)   
            
            dicecoeff = (2 * area_intersect) / (area_ct + area_pet)
            dice.append(dicecoeff)

            z.append(j)    
            

        # Get ref file for calculate distance between CT and PET.
        Ref_CT=nrrd.read(os.path.join(pathtemp, 'VOI_' + modality + '_' + names + '_resampled.nrrd'))[1]
        Ref_PET=nrrd.read(os.path.join(pathtemp, 'VOI_PET_' + names + '_con_cabeceraPET.nrrd'))[1]
        

        # Load position left corner pixel
        PixelPosition_CT = (float(Ref_CT['space origin'][0]), float(Ref_CT['space origin'][1]), float(Ref_CT['space origin'][2]))
        PixelPosition_PET = (float(Ref_PET['space origin'][0]), float(Ref_PET['space origin'][1]), float(Ref_PET['space origin'][2]))
        ConstPixelSpacing_CT = (float(Ref_CT['space directions'][0][0]), float(Ref_CT['space directions'][1][1]), float(Ref_CT['space directions'][2][2]))
        ConstPixelSpacing_PET = (float(Ref_PET['space directions'][0][0]), float(Ref_PET['space directions'][1][1]), float(Ref_PET['space directions'][2][2]))


        # Calculo en mm del primer corte del CT y del PET y su diferencia.
        pos_PET=PixelPosition_PET[2]+ConstPixelSpacing_PET[2]*li_pet
        pos_CT=PixelPosition_CT[2]+ConstPixelSpacing_CT[2]*li_ct
        DIF_CORTES=abs(pos_PET-pos_CT)
 
                       
    else:

        if phantom == 'JASZCZAK':
            if seg_desired == '1':
                i=0
            elif seg_desired == '2':
                i=1
            elif seg_desired == '3':
                i=2
            elif seg_desired == '4':
                i=3
            elif seg_desired == '5':
                i=4
            elif seg_desired == '6':
                i=5

        if phantom == '3D':
            if seg_desired == 'right':
                i = 2
            elif seg_desired == 'center':
                i = 1
            elif seg_desired == 'left':
                i = 0



            
        nrrd.write('Seg_ref_' + names[i] + '.nrrd', VOIs_REF[i], header_PET, index_order='F')
        nrrd.write('VOI_PET_' + names[i] + '.nrrd', VOIs_PET[i], index_order='F')
        nrrd.write('VOI_PET_' + names[i] + '_con_cabeceraPET.nrrd', VOIs_PET, header_PET, index_order='F')
        VOI_PET_cabecera = nrrd.read(os.path.join(pathtemp, 'VOI_PET_' + names[i] + '_con_cabeceraPET.nrrd'))[0]



        # Intersection TH+VOI
        ARR_functions.resampling(path_CT, os.path.join(pathtemp, 'Seg_ref_' + names[i] + '.nrrd'), pathtemp, 'VOI_REF_' + modality + '_' + names[i])        # VOI CT resampling

        
        ref_CT = nrrd.read(os.path.join(pathtemp, 'VOI_REF_' + modality + '_' + names[i] + '_resampled.nrrd'))[0].astype('bool')
        th_CT = np.copy(Threshold_CT)


        th_CT[~ref_CT] = 0
        th_CT.astype('uint8')
        nrrd.write('VOI_CT_' + names[i] + '.nrrd', th_CT, header_CT, index_order='F')

        # Intersection TH+VOI CT
        ARR_functions.resampling(path_PET, os.path.join(pathtemp, 'VOI_' + modality + '_' + names[i] + '.nrrd'), pathtemp,
                                 'VOI_' + modality + '_' + names[i])                                                                               # VOI CT resampling
        
        # Intersection VOI_PET+VOI_CT
        VOI_CT = nrrd.read(os.path.join(pathtemp, 'VOI_' + modality + '_' + names[i] + '_resampled.nrrd'))[0].astype('bool')
        intersect = VOI_CT * VOI_PET_cabecera
        nrrd.write('intersection_' + names[i] + '.nrrd', intersect, header_PET, index_order='F')




        # Number of slices calculation
        li_pet= 0
        li_ct=0

       
        for s in range(0, len(VOI_PET_cabecera[i][0, 0, :])):
            if np.mean(VOI_PET_cabecera[i][:, :, s]) == 0:
                li_pet = li_pet+1
            else:
                break
        
        for s in range(li_pet, len(VOI_PET_cabecera[i][0, 0, :])):
            if np.mean(VOI_PET_cabecera[i][:, :, s]) != 0 :      
                lf_pet = s
            else:
                break
        print('La imagen PET va desde el corte ', li_pet, 'hasta ', lf_pet)
        
       
        

        for r in range(0, len(VOI_CT[0, 0, :])):
            if np.mean(VOI_CT[:, :, r]) == 0:
                li_ct = li_ct + 1
            else:
                break


        for s in range(li_ct, len(VOI_CT[0, 0, :])):
            if np.mean(VOI_CT[:, :, s])!= 0 :
                lf_ct = s 
            else:
                break
        print('La imagen CT va desde el corte ', li_ct, 'hasta ', lf_ct)




        
        # Elección de los cortes inicial y final para calcular Dice
       
        if li_pet<li_ct:    
            li=li_pet
            inc_z=li_ct-li_pet
        else:               
            li=li_ct
            inc_z=li_pet-li_ct

        if lf_pet>lf_ct:
            lf=lf_pet
        else:
            lf=lf_ct


        
        # DICE calculation
        dice = []
        z = []

        for j in range(li, lf+1):       
            ct_slice = VOI_CT[:, :, j]
            pet_slice = VOI_PET_cabecera[i][:, :, j]
            intersect_slice = intersect[:, :, j]
            area_ct = np.sum(ct_slice == True)
            area_pet = np.sum(pet_slice == 1)   
            area_intersect = np.sum(intersect_slice == 1)   
            
            dicecoeff = (2 * area_intersect) / (area_ct + area_pet)
            dice.append(dicecoeff)

            z.append(j)    


        # Get ref file for calculate distance between CT and PET.
        Ref_CT=nrrd.read(os.path.join(pathtemp, 'VOI_' + modality + '_' + names[i] + '_resampled.nrrd'))[1]
        Ref_PET=nrrd.read(os.path.join(pathtemp, 'VOI_PET_' + names[i] + '_con_cabeceraPET.nrrd'))[1]
        

        # Load position left corner pixel
        PixelPosition_CT = (float(Ref_CT['space origin'][0]), float(Ref_CT['space origin'][1]), float(Ref_CT['space origin'][2]))
        PixelPosition_PET = (float(Ref_PET['space origin'][0]), float(Ref_PET['space origin'][1]), float(Ref_PET['space origin'][2]))
        ConstPixelSpacing_CT = (float(Ref_CT['space directions'][0][0]), float(Ref_CT['space directions'][1][1]), float(Ref_CT['space directions'][2][2]))
        ConstPixelSpacing_PET = (float(Ref_PET['space directions'][0][0]), float(Ref_PET['space directions'][1][1]), float(Ref_PET['space directions'][2][2]))


        # Calculo en mm del primer corte del CT y del PET y su diferencia.
        pos_PET=PixelPosition_PET[2]+ConstPixelSpacing_PET[2]*li_pet
        pos_CT=PixelPosition_CT[2]+ConstPixelSpacing_CT[2]*li_ct
        DIF_CORTES=abs(pos_PET-pos_CT)





    mean_dice = statistics.mean(dice)
    StandarDesviation_dice=statistics.stdev(dice)
    print (mean_dice, '\u00B1',  StandarDesviation_dice)
    

    #print (dice)   
 
    # DICE's plot according to the slice
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(111)

    ax1.plot(z, dice, color=cls.to_rgba('C7', 0.7))
    ax1.axhline(y=mean_dice, color=cls.to_rgba('C9', 0.8))
    
    if phantom == 'JASZCZAK':
        ylimit=0.68             # limit value for coregistration
     
    else:
        ylimit=0.93             # limit value for coregistration

    ax1.axhline(y=ylimit, color=cls.to_rgba('#FFA500', 0.8), linestyle='dotted')       
    ax1.fill_between(z, dice, color=cls.to_rgba('C7', 0.3))

    ax1.set_xlim([min(z), max(z)])
    ax1.set_ylim([0, 1])    
    ax1.set_xlabel('Slices', fontsize=15)
    ax1.set_ylabel('DICE', fontsize=15)
    if activity == 'yes':
        ax1.set_title('Dice similarity index ' + names[i] + ' tube', fontsize=20)
    else:
        ax1.set_title('Dice similarity index ', fontsize=20)
    
    if round(mean_dice, 2) >= ylimit:   
        if DIF_CORTES>2:
            color=('black','red', 'black')
            ax1.legend([f'DICE coefficient',f'Mean DICE= {round(mean_dice, 2)} \u00B1 {round(StandarDesviation_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm \nThe PET ' + modality + ' images are not well co-registered',
                    f'DICE threshold value: ' + str(ylimit)],labelcolor=color,title=f'Co-registration properties')
        else:
            color=('black','green', 'black')
            ax1.legend([f'DICE coefficient', f'Mean DICE= {round(mean_dice, 2)} \u00B1 {round(StandarDesviation_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm \nThe PET ' + modality + ' images are well co-registered',
                    f'DICE threshold value: ' + str(ylimit)], labelcolor=color, title=f'Co-registration properties')
    else:
        color=('black','red', 'black')
        ax1.legend([f'DICE coefficient']+ [f'Mean DICE= {round(mean_dice, 2)} \u00B1 {round(StandarDesviation_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm \nThe PET ' + modality + ' images are not properly co-registered',
             f'DICE threshold value: ' + str(ylimit)], labelcolor=color, title=f'Co-registration properties')        


    pdf.savefig(fig)
    pdf.close()
    print ('fin del programa')

