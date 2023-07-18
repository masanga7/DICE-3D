import math
import os
import numpy as np
import nrrd
import cc3d
import re
import pydicom as dicom
import pandas as pd
import matplotlib.pyplot as plt #linea nueva
import pydicom                  #linea nueva
from natsort import natsorted
from os import walk
from os.path import splitext, join
from skimage import measure, morphology
from skimage.feature import canny
from scipy import ndimage
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog




def thresholdseg(image, th_value, phantom, size_voxel):                                                             
	imageseg = np.copy(image)                                                                       
	if th_value == 'water CT':                                                                      
                if phantom == '3D':
                        th = np.where((imageseg >= -50) & (imageseg <= 500), 1, 0)                      
                        fill = ndimage.binary_fill_holes(th).astype(int)                                
                        erosion = morphology.binary_erosion(fill)                                       
                        cleared = morphology.binary_erosion(erosion).astype(np.uint8)                   
                
                        
                if phantom == 'JASZCZAK':

                        print ('has entrado en CT con el maniquí Jacksack')
                        vmax = np.nanmax(imageseg)                  
                        vmin=np.nanmin(imageseg)
                        vtotal=vmax+abs(vmin)
                        value_th = 0.55
                        th = np.where((imageseg >= value_th*vtotal - abs(vmin)), 1, 0)

                        for s in range (0,10):                                                          #Para poner la base del maniquí en negro
                                th[:,:,s] = 0

                        ############ Para saber el centro del maniquí#############
                        s=0
                        s_z=int(len(imageseg[0,0,:])/2)
                        s_x = int(len (imageseg[:,0,0])/2)
                        s_y = 0
                        size_pixel = float(size_voxel[0][0])
                        print ('size_pixel', size_pixel)
                        slice_thickness = float(size_voxel[2][2])
                        print ('thickness', slice_thickness)
                        dim_phan=int(192.2/(2*size_pixel))                                              #192.2 es el diametro del maniquí
                        
                        for s_y in range (0, len(imageseg[s_x,:,s_z])):
                                if imageseg[s_x,s_y, s_z]>=0:
                                        print ('El borde del maniquí se encuentra en. ', s_y)
                                        print ('El centro del maniquí se encuentra en: ', s_x, s_y + dim_phan, s_z)
                                        break

                        ########### Para eliminar los voxeles blancos que están a una distancia mayor que Rmaniqui-d8pixeles  ####
                        z=0
                        y=0
                        x=0
                        for z in range (0, len(th[0,0,:])):
                                for y in range (0, len(th[:,:,:])):
                                        for x in range (0, len (th[:,:,:])):
                                                dif_pos=((1.171875*(x-s_x))**2+(1.171875*(y-(s_y + dim_phan)))**2)**0.5                 
                                                if (dif_pos>=86.725):                                                                   
                                                        th[x,y,z]=0




                        ######### Area dentro de las esferas ##########
                        
                        fill = ndimage.binary_fill_holes(th).astype(int)                                                
                        mis_esferas1 = np.where((fill==1) & (th==1),0,fill)                                                                      
                        


                        ######### Reduciendo esfera #########

                        value_th1 = 0.58
                        th1=np.where(imageseg >= value_th1*vtotal - abs(vmin),1,0)                                                        
                        print ('Aplicando después al CT un porcentaje del ', value_th1*100, 'es: ', value_th1*vtotal - abs(vmin))

                        intersec = np.where((fill==1) & (th1==1),0,fill)                                                
                        erosion= morphology.binary_erosion(intersec)                        
                        mis_esferas=ndimage.binary_dilation(erosion).astype(int)                                        
                        fill1 = np.where((mis_esferas==1) & (th1==0),1,th1)                                                                     


                        ######## Creando Volumen de las esferas para ordenarlas en el sentido de las agujas del reloj (falta terminar)#######
                        labels_in, N = cc3d.connected_components(mis_esferas, connectivity=26, return_N=True)           

                        print ('Número de muestras detectadas',N)
                        images = []                                                                                     
                        centroids_CT = []                                                                               
                        voxels = []
                        for segid in range(1, N + 1):                                                                   
                                extracted_image = labels_in * (labels_in == segid)                                      
                                extracted_image[extracted_image != 0] = 1                                               
                                stats = cc3d.statistics(extracted_image)                                               
                                centroids_CT.append(stats['centroids'][-1])                                             
                                voxels.append(stats['voxel_counts'][-1])
                                print ('numero de voxels', voxels[segid-1])
                                images.append(list(extracted_image))
                                print (segid, ' centroide_CT: ', centroids_CT[segid-1])

                        pairs = list(zip(images, centroids_CT, voxels))                                                 
                        pairs.sort(key=lambda x: x[1][1])                                                               
                        images, centroids_CT, voxels = map(np.array, zip(*pairs))                                       

                        cleared=mis_esferas
                        

	if th_value == 'water MRT2':
                print ('has entrado en MRT2')
                vmax = np.nanmax(imageseg)                  
                vmin=np.nanmin(imageseg)
                vtotal=vmax+abs(vmin)
                value_th = 0.30
                th = np.where((imageseg >= value_th*vtotal - abs(vmin)), 1, 0)
                cleared=th


                
	elif th_value == 'water MRT1':
                print ('has entrado en MRT1')

                vmax = np.nanmax(imageseg)                  
                vmin=np.nanmin(imageseg)
                vtotal=vmax+abs(vmin)
                value_th = 0.30
                print ('El valor absoluto minimo aplicando a la RM un porcentaje del ', value_th*100, 'es: ', value_th*vtotal - abs(vmin))

                th = np.where((imageseg >= value_th*vtotal - abs(vmin)), 1, 0)                                                                  
                fill = ndimage.binary_fill_holes(th, np.ones((5,5,1))).astype(int)
                erosion = morphology.binary_erosion(fill, np.ones((2,2,1))).astype(int)
                dilated=ndimage.binary_dilation(erosion, np.ones((2,2,1))).astype(int)

                labels_in, N = cc3d.connected_components(dilated, connectivity=26, return_N=True)                                               
                print ('Número de muestras detectadas',N)
                cleared =labels_in
                
                
	return cleared                                                                  




def reference_seg(image, threshold, phantom):                                                                                           

	imageseg = np.copy(image)                                                                                                       
	vmax = np.nanmax(imageseg)                                                                                                      
	vmin=np.nanmin(imageseg)
	vtotal=vmax+abs(vmin)
	if phantom == 'JASZCZAK':
                for z in range (int(len (imageseg[0,0,:]))-5, len (imageseg[0,0,:])):
                        for y in range (0, len (imageseg[:,:,:])):
                                for x in range (0, len (imageseg[:,:,:])):
                                        imageseg[x,y,z] = 0
                hole_coords=np.where(imageseg[:, :, :] >= threshold*vtotal - abs(vmin))                                              
	else:
                hole_coords = np.where(imageseg[:, :360, :] >= threshold*vtotal - abs(vmin))                                             
	
	imageseg2 = np.copy(image) * 0                                                                                                  
	imageseg2[hole_coords] = 1                                                                                                      
	erosion = morphology.binary_erosion(imageseg2)                                                                                  
	dilated = ndimage.binary_dilation(erosion)                                                                                      
	dilated = dilated.astype(np.int16)
	
	labels_in, N = cc3d.connected_components(dilated, connectivity=26, return_N=True)                                               
	print ('Número de muestras detectadas',N)
	images = []                                                                                                                                                                  
	images_ref = []                                                                                                                 
	centroids = []                                                                                                                  
	for segid in range(1, N + 1):                                                                                                   
		extracted_image = labels_in * (labels_in == segid)                                                                      
		extracted_image[extracted_image != 0] = 1                                                                               
		stats = cc3d.statistics(extracted_image)                                                                                

		if phantom == 'JASZCZAK':
                        centroids.append(stats['centroids'][-1])                                                                        
                        extracted_image_dilated = ndimage.binary_dilation(extracted_image, iterations=8).astype(np.uint16)
		else:
                        centroids.append(stats['centroids'][-1][0])                                                                     
                        extracted_image_dilated = ndimage.binary_dilation(extracted_image, iterations=8).astype(np.uint8) 
		images.append(list(extracted_image))
		images_ref.append(list(extracted_image_dilated))
		print (segid, ' centroide: ', centroids[segid-1])

	if phantom == 'JASZCZAK':
                pairs = list(zip(images, images_ref, centroids))                                                                       
                pairs.sort(key=lambda x: x[2][1])                                                                                       
                images, images_ref, centroids = map(np.array, zip(*pairs))                                                              

	return images, images_ref, centroids                                                                                            




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
