#search folders for vis arm files to be telluric-corrected
#place inside end products folder for a given target, should have VIS subfolder
import os

if (not 'TELCORR' in os.listdir()):
    os.system('mkdir TELCORR')

folders_list = [f for f in os.listdir('VIS/') if f != 'README']
folders_list.sort()

for i in range(len(folders_list)):
    index = '_'+str(i+1).zfill(2)
    path = 'VIS/'+folders_list[i]
    file = [f for f in os.listdir(path) if '_IDP_VIS.fits' in f][0]
    new_filename = file[:-5]+index+'.fits'
    os.system(f'cp {path+'/'+file} TELCORR/{new_filename}')