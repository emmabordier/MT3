import numpy as np
import random
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits
import os as os
from astropy import wcs
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from scipy.optimize import curve_fit

f=[100,143,217,353,545,857]

def writefits(array, fname, overwrite=False):

  if (os.path.isfile(fname)):
    if (overwrite == False):
      print("File already exists!, skiping")
      return
    else:
      print("File already exists!, overwriting")
      os.remove(fname)
  try: 
    hdu = fits.PrimaryHDU(array)
    hdu.writeto(fname)
  except:
    print('b1') 
    hdu = pf.PrimaryHDU(fname)
    print('b2') 
    hdu.writeto(array)

  print("File: %s saved!" % fname)

  return



def reduceimage(f): #f=frequence
    NSIDE=2048
    map=hp.read_map('/home/vbonjean/PLANCK/HFI_SkyMap_'+str(f)+'_2048_R3.01_full.fits')
    mapred=hp.ud_grade(map,256)
    hp.write_map('HFI_SkyMap_reduced_'+str(f)+'.fits', mapred)

def plotimage(f):
  image=hp.read_map('../data/HFI_SkyMap_reduced_'+str(f)+'.fits')
  hp.mollview(image, min=0, max=1000)
  plt.show()


#construction patch


def myWCS(N,lon,lat,tp): #tp taille pixel 1.7 lon lat coord point central
    w=wcs.WCS(naxis=2)
    w.wcs.crpix=[N/2,N/2]
    w.wcs.cdelt=[-tp/60,tp/60]
    w.wcs.crval=[lon,lat]
    w.wcs.ctype=["GLON-TAN","GLAT-TAN"]
    return(w)
      

def patch(planckmap,N,lon,lat,tp,nside): #N=nombre de pixels patch tp=taille pix arcmin
    w=myWCS(N,lon,lat,tp)
    patch=np.zeros((N,N))
    i,j=np.indices((N,N))
    coord=w.wcs_pix2world(i,j,1)
    ipix=hp.ang2pix(nside,coord[0],coord[1], lonlat=True) #nest=True lonlat=True pour degres
    patch=planckmap[ipix]
    return(patch)
            

"""
A=patch(100,300,0,0,1.7,256)
B=patch(100,40,0,0,14.2,256) #cf calcul

plt.figure()
plt.subplot(1,2,1)
plt.title('Test image apres transformation pixels \n (taille pixel=1.7arcmin)')
plt.imshow(A)
plt.subplot(1,2,2)
plt.title('Image apres transformation pixels \n  (taille du pixel =14.2arcmin)')
plt.imshow(B)
plt.show()

"""

#Ouverture du catalogue et extraction des coordonnees des sources
Catalogue=fits.open('../data/HFI_PCCS_SZ-union_R2.08.fits')
hdu=Catalogue[1]
data=hdu.data
Glon=data['GLON']
Glat=data['GLAT']

'''
def plotpatch(N,lon,lat,tp,nside,S): #S est la source choisie
    M=[]
    f=[100,143,217,353,545,857]
    for i in f:
        M.append(patch(i,40,Glon[S],Glat[S],14.2,256))
    return(M)
'''

def unites_SZ(f):
    L1=[100,143,217,353]
    L2=[545,857]
    ucsz=[-0.24815,-0.35923,5.152,0.161098,0.06918,0.0380]
    uciras=[244.1,371.74,483.690,287.450,58.04,2.27]
    if f in L1:
        i=L1.index(f)
        map=hp.read_map('/home/vbonjean/PLANCK/HFI_SkyMap_'+str(f)+'_2048_R3.01_full.fits')
        map=ucsz[i]*map
        hp.write_map('HFI_calibre_'+str(f)+'.fits', map)
    if f in L2:
        i=L2.index(f)+4
        map=hp.read_map('/home/vbonjean/PLANCK/HFI_SkyMap_'+str(f)+'_2048_R3.01_full.fits')
        map=ucsz[i]/uciras[i]*map
        hp.write_map('HFI_calibre_'+str(f)+'.fits', map)
    if f==70:
      map=hp.read_map('../data/LFI_SkyMap_0'+str(f)+'-BPassCorrected-field-IQU_1024_R3.00_full.fits')
      map=ucsz70*map
      hp.write_map('LFI_calibre_'+str(f)+'.fits', map)

def unites_KCMB(f):
    L1=[100,143,217,353]
    L2=[545,857]
    uciras=[58.04,2.27]
    if f in L1:
        map=hp.read_map('/home/vbonjean/PLANCK/HFI_SkyMap_'+str(f)+'_2048_R3.01_full.fits')
        hp.write_map('HFI_calibre_KCMB_'+str(f)+'.fits', map)
    if f in L2:
        i=L2.index(f)
        map=hp.read_map('/home/vbonjean/PLANCK/HFI_SkyMap_'+str(f)+'_2048_R3.01_full.fits')
        map=map/uciras[i]
        hp.write_map('HFI_calibre_KCMB_'+str(f)+'.fits', map)

def unites_MJY(f):
    L1=[100,143,217,353]
    L2=[545,857]
    uciras=[244.1,371.74,483.69,287.45]
    if f in L1:
        i=L1.index(f)
        map=hp.read_map('/home/vbonjean/PLANCK/HFI_SkyMap_'+str(f)+'_2048_R3.01_full.fits')
        map=map*uciras[i]
        hp.write_map('HFI_calibre_MJY_'+str(f)+'.fits', map)
    if f in L2:
        map=hp.read_map('/home/vbonjean/PLANCK/HFI_SkyMap_'+str(f)+'_2048_R3.01_full.fits')
        hp.write_map('HFI_calibre_MJY_'+str(f)+'.fits', map)

  
"""
#changer unites

for i in f:
  unites_KCMB(i) #unites_SZ(i)  unites_MJY(i)

"""
  
"""      
#plot patch
  
M=plotpatch(40,Glon[100],Glat[100],14.2,256,100)

plt.figure()
plt.suptitle('Test patch source 100 pour les differentes frequences')
plt.subplot(2,3,1)
plt.imshow(M[0])
plt.title('f=100')
plt.subplot(2,3,2)
plt.imshow(M[1])
plt.title('f=143')
plt.subplot(2,3,3)
plt.imshow(M[2])
plt.title('f=217')
plt.subplot(2,3,4)
plt.imshow(M[3])
plt.title('f=353')
plt.subplot(2,3,5)
plt.imshow(M[4])
plt.title('f=545')
plt.subplot(2,3,6)
plt.imshow(M[5])
plt.title('f=857')
plt.show()
"""

#ecrire fichiers patchs

Catalogue2=np.loadtxt('../data/potential_amas_ana_emma.txt')
longitude=Catalogue2[:,0]
latitude=Catalogue2[:,1]

def freqpatch(f,N,tp,nside,longitude, latitude, nom): #N:nombre de pixels patch longitude latitude du catalogue choisi (pccs ou Hubert Stephano)
  image=hp.read_map('../data/HFI_calibre_MJY_'+str(f)+'.fits') #changer
  Tableau=np.zeros((len(longitude),N,N))
  for i in range(len(longitude)):
    Tableau[i,:,:]=patch(image,N,longitude[i],latitude[i],tp,nside)
  writefits(Tableau, 'Patch_'+nom+'_'+str(f)+'.fits') #changer nom

def patch_confirmed(f):
  data=fits.open('../data/HFI_PCCS_SZ-union_R2.08.fits')[1].data
  r=data['REDSHIFT']
  amas=np.where(r>0)
  carte=fits.open('../data/Patch_MJY_'+str(f)+'.fits')[0].data #changer 
  carte_conf=[]
  for i in amas[0]:
    carte_conf.append(carte[i,:,:])
  carte_conf=np.array(carte_conf)
  writefits(carte_conf, 'Patch_MJY_confirmed_amas_'+str(f)+'.fits') #changer nom

def patch_candidats(f):
  data=fits.open('../data/HFI_PCCS_SZ-union_R2.08.fits')[1].data
  q=data['Q_NEURAL']
  candidats=np.where(q>0.7)  #il y en a 1447 sources
  carte=fits.open('../data/Patch_MJY_'+str(f)+'.fits')[0].data #changer 
  carte_cand=[]
  for i in candidats[0]:
    carte_cand.append(carte[i,:,:])
  carte_cand=np.array(carte_cand)
  writefits(carte_cand, 'Patch_MJY_candidats_amas_'+str(f)+'.fits') #changer nom

"""
for i in f:
  freqpatch(i,128,1.5,2048)

for i in f:
  patch_confirmed(i)

for i in f:
  patch_candidats(f)

"""

"""    
#comparer patch a gnomview

carte=hp.read_map('../data/HFI_calibre_100.fits')

test=patch2(carte,256,Glon[10],Glat[10],1.7,2048)

hp.gnomview(carte,rot=[Glon[10],Glat[10]],xsize=256,ysize=256) #changer a 128
plt.figure()
plt.imshow(test.T, origin='lower')
plt.show()
"""


def plotpatchfreq(nom,S): #S numero de la source voulue
  f=[100,143,217,353,545,857]
  L=[]
  for i in f:
    carte=fits.open('../data/Patch_'+nom+'_'+str(i)+'.fits')[0].data #changer pour SZ
    L.append(carte[:,:,S])
  plt.figure()
  plt.suptitle('Patchs source'+str(S)+' du catalogue a différentes fréquences')
  plt.subplot(2,3,1)
  plt.imshow(L[0].T,origin='lower')
  plt.title('Patch pour f=100GHz')
  plt.subplot(2,3,2)
  plt.imshow(L[1].T,origin='lower')
  plt.title('Patch pour f=143GHz')
  plt.subplot(2,3,3)
  plt.imshow(L[2].T,origin='lower')
  plt.title('Patch pour f=217GHz')
  plt.subplot(2,3,4)
  plt.imshow(L[3].T,origin='lower')
  plt.title('Patch pour f=353GHz')
  plt.subplot(2,3,5)
  plt.imshow(L[4].T,origin='lower')
  plt.title('Patch pour f=545GHz')
  plt.subplot(2,3,6)
  plt.imshow(L[5].T,origin='lower')
  plt.title('Patch pour f=857GHz')
  plt.tight_layout()
  plt.show()


def stack(nom,f):
  carte=fits.open('../data/Patch_'+nom+'_'+str(f)+'.fits')[0].data #changer 
  stack=np.mean(carte,axis=0)
  writefits(stack, 'Stack_'+nom+'_'+str(f)+'.fits') #changer

def stack_confirmed(f):
  data=fits.open('../data/HFI_PCCS_SZ-union_R2.08.fits')[1].data
  r=data['REDSHIFT']
  amas=np.where(r>0)
  #q=data['Q_NEURAL']
  #index=np.where(q>0.95)
  carte=fits.open('../data/Patch_KCMB_'+str(f)+'.fits')[0].data #changer 
  carte_conf=[]
  for i in amas[0]:
    carte_conf.append(carte[i,:,:])
  carte_conf=np.array(carte_conf)
  stack=np.mean(carte_conf,axis=0)
  #return np.shape(carte_conf)
  writefits(stack, 'Stack_KCMB_confirmed_amas_'+str(f)+'.fits') #changer nom

def stack_bins(f): #pour 10 bins
  data=fits.open('../data/HFI_PCCS_SZ-union_R2.08.fits')[1].data
  q=data['Q_NEURAL']
  qord=np.argsort(q)
  carte=fits.open('../data/Patch_MJY_'+str(f)+'.fits')[0].data
  carte_bins=[]
  stack=[]
  for i in qord:
    carte_bins.append(carte[i,:,:])
  carte_bins=np.array(carte_bins)
  for j in range(0,1650,165):
    coupe=carte_bins[j:j+165,:,:]
    stack.append(np.mean(coupe, axis=0))
  writefits(stack, 'Stack_MJY_bins_q_'+str(f)+'.fits')

"""
for i in f:
  stack(i)

for i in f:
  stack_confirmed(i)
""" 

""" 
#Pixel chaud

s=fits.open('../data/Stack_100.fits')[0].data  
a=np.where(s==np.max(s))
i=a[0][0]
j=a[1][0]

s[i,j]=0
"""

def plotstacks(nom): 
  f=[100,143,217,353,545,857]
  L=[]
  for i in f:
    s=fits.open('../data/Stack_'+nom+'_'+str(i)+'.fits')[0].data #changer
    L.append(s)
  plt.figure()
  plt.suptitle('Stacks 44 candidats Hubert et Stephano') #1094 amas confirmes 1447 pour q>0.7
  plt.subplot(2,3,1)
  plt.imshow(L[0].T,origin='lower')
  plt.title('Stack pour f=100GHz')
  plt.subplot(2,3,2)
  plt.imshow(L[1].T,origin='lower')
  plt.title('Stack pour f=143GHz')
  plt.subplot(2,3,3)
  plt.imshow(L[2].T,origin='lower')
  plt.title('Stack pour f=217GHz')
  plt.subplot(2,3,4)
  plt.imshow(L[3].T,origin='lower')
  plt.title('Stack pour f=353GHz')
  plt.subplot(2,3,5)
  plt.imshow(L[4].T,origin='lower')
  plt.title('Stack pour f=545GHz')
  plt.subplot(2,3,6)
  plt.imshow(L[5].T,origin='lower')
  plt.title('Stack pour f=857GHz')
  plt.tight_layout()
  plt.show()



"""
#largeur a mi-hauteur "source"

def gauss(x,a,xo,sigma):
  return a*np.exp(-(x-xo)**2/(2*sigma**2))


def fwhm(f):
  stack=fits.open('../data/Stack_128_'+str(f)+'.fits')[0].data
  coupe=stack[64,:]
  coupe=coupe-np.min(coupe)
  sscoupe=coupe[40:88]
  sscoupe=sscoupe-np.min(sscoupe)
  x=np.arange(128)
  x2=np.arange(40,88,1)
  a=np.max(coupe)-np.min(coupe)
  popt1,pcov1=curve_fit(gauss,x2,sscoupe,p0=(a,24,10))
  fwhm=popt1[2]*np.sqrt(2*np.log(2))*2
  plt.figure()
  plt.plot(x2,sscoupe)
  plt.plot(x2,gauss(x2,*popt1),label='fit: fwhm='+str(fwhm))
  plt.legend()
  plt.title('Largeur a mi-hauteur stack source pour f='+str(f)+'GHz')
  plt.show()
"""
  
def astrobootstrap(Nb,f):
  data=fits.open('../data/Patch_MJY_confirmed_amas_'+str(f)+'.fits')[0].data
  with NumpyRNGContext(1):
    bootresult=bootstrap(data,Nb) #array de dim [Nb,S,N,N]
  writefits(bootresult, 'bootstrap_MJY_'+str(f)+'.fits')

def stackalea(f,nom): #nom=MJY_candidats_amas ou MJY_confirmed_amas
  data=fits.open('../data/Patch_'+nom+'_'+str(f)+'.fits')[0].data 
  L=np.arange(len(data))
  Stack=np.zeros((128,128)) #stack avec echantillons aleatoires
  for i in range(len(data)):
    indice=random.choice(L)
    Stack+=data[indice,:,:]
  return Stack 
    
def mybootstrap(Nb,f,nom):
  boot=[]
  for i in range(Nb):
    boot.append(stackalea(f,nom))
  boot=np.array(boot)
  boot=boot/Nb
  writefits(boot, 'bootstrap_'+nom+'_'+str(f)+'_'+str(Nb)+'.fits')
    
  
    
    
