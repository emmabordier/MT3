import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits
import os as os
from astropy import wcs
from scipy.optimize import curve_fit

freq=np.array([100,143,217,353,545,857])

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


def round_mask(shape, radius, x_0=None, y_0=None):
    """ Generate a a round mask centered in x_0, y_0.

    Parameters
    ==========
    shape : tuple
        A 2-elements tuple containing the shape of the mask. 
    radius : int
        The radius of the mask

    Returns
    =======
    A np.ndarray of `shape`, where all pixels that are closer than `radius` to
    the center are set at 1, and all other pixels to 0. 

    Example
    =======
    >>> round_mask((4,4), 2)
    array([[ 0.,  0.,  0.,  0.],
           [ 0.,  0.,  1.,  0.],
           [ 0.,  1.,  1.,  1.],
           [ 0.,  0.,  1.,  0.]])
    """
    if not x_0:
        x_0 = shape[0]//2
    if not y_0:
        y_0 = shape[1]//2 

    mask = np.ndarray(shape)
    i_0 = x_0
    j_0 = y_0
    for i in range(shape[0]):
        for j in range(shape[1]):
            if np.sqrt((i-i_0)**2 + (j-j_0)**2) <= radius:
                mask[i,j] = 1
            else:
                mask[i,j] = 0
    return mask


def aperture_photom(data, peaks, fwhm, fr1,fr2,fr3):
    """ Compute aperture photometry on data 

    Parameters
    ==========
    data : 
        a 2D np.ndarray containing the *original* image
    peaks : 
        a 2D np.ndarray of booleans marking the position of sources, same dimensions as 'data'
    fwhm : 
        the fwhm (in px) of the PSF

    Return a tuple containing:
        - source (x, y) (in px)
        - intensity (in DN)
        - noise 
        - SN
    """

    frame_size = 14*int(fwhm)

    # Radius of the disk
    r1 = fr1 * fwhm  #r1 = 1.5 * fwhm 
    disk = round_mask((frame_size, frame_size), r1)

    # Inner and outer radii of the annulus
    r2 = fr2 * fwhm #r2 = 2.5 * fwhm
    r3 = fr3 * fwhm #r3 = 3 * fwhm
    corona = round_mask((frame_size, frame_size), r3)-round_mask((frame_size, frame_size), r2)

    surface_ratio = disk.sum() / corona.sum()
    sources_catalog = []

    for i, j in np.array(np.where(peaks)).T:
        tranche = data[i-frame_size//2:i+frame_size//2, j-frame_size//2:j+frame_size//2]
        flux = np.sum(disk * tranche)
        background = surface_ratio * np.sum(corona * tranche)
        noise = surface_ratio * np.std(corona * tranche)
        intensity = (flux - background)
        SN = flux / noise
        sources_catalog.append(((i, j), intensity, noise, SN))
    return (sources_catalog,surface_ratio,disk,corona,tranche)


"""
        except ValueError:
            # This happens when either the disk or annulus contains at least a
            # single nan
            try:
                mask_c = np.isnan(corona * tranche)
                mask_d = np.isnan(disk * tranche)
                if (mask_c.sum() <= 10 and mask_d.sum() <= 2):
                    corona = np.ma.array(corona, mask=mask_c).filled(fill_value=0)
                    disk = np.ma.array(disk, mask=mask_d).filled(fill_value=0)
                    surface_ratio = disk.sum() / corona.sum()
                    flux = np.sum(disk * tranche)
                    background = surface_ratio * np.sum(corona * tranche)
                    noise = surface_ratio * np.std(corona * tranche)
                    intensity = (flux - background)
                    SN = flux / noise
                    sources_catalog.append(((i, j), intensity, noise, SN))
                else:
                    raise ValueError
            except ValueError:
                sources_catalog.append(((i, j), np.nan, np.nan, np.nan))
"""
   
def photomstack(f,fwhm,nom, fr1, fr2, fr3): #2.5 3 4 (5,6)
  stack=fits.open('../data/Stack_'+nom+'_'+str(f)+'.fits')[0].data # changer 
  shape=np.shape(stack)
  peaks = np.zeros(shape, dtype=bool)
  peaks[64,64]=True
  #photom = aperture_photom(stack, peaks, fwhm, 1.5, 2.5, 3) #choisir valeurs fr
  photom = aperture_photom(stack, peaks, fwhm, fr1, fr2, fr3) # 5 6
  return photom


#photometrie pour les 6 stacks

def psf(tp): #taille pixel
  FWHMa=np.array([9.66,7.22,4.90,4.92,4.67,4.22]) #en arcmin  #13.21 pour f=70
  FWHMp=FWHMa/tp #en pixels
  return FWHMp

RES_PHOTOM4=[]
photom4=[] #3 indices avec differentes photometries
noise4=[]
RES_PHOTOM5=[]
photom5=[] 
noise5=[]
RES_PHOTOM6=[]
photom6=[] 
noise6=[]
for i in range(6):
  f=freq[i]
  fwhm=psf(1.5)[i] #1.5= taille du pixel
  RES_PHOTOM4.append(photomstack(f,fwhm,'MJY_confirmed_amas',2.5,3,4)) 
  photom4.append(RES_PHOTOM4[i][0][0][1])
  noise4.append(RES_PHOTOM4[i][0][0][2])
  RES_PHOTOM5.append(photomstack(f,fwhm,'MJY_confirmed_amas',2.5,3,5)) 
  photom5.append(RES_PHOTOM5[i][0][0][1])
  noise5.append(RES_PHOTOM5[i][0][0][2])
  RES_PHOTOM6.append(photomstack(f,fwhm,'MJY_confirmed_amas',2.5,3,6)) 
  photom6.append(RES_PHOTOM6[i][0][0][1])
  noise6.append(RES_PHOTOM6[i][0][0][2])


#photometrie candidats q>0.7

def photometrie(nom): #fonction retournant photometrie pour differentes frequences
  photom=[]
  for i in range(6):
    f=freq[i]
    fwhm=psf(1.5)[i]
    PHOTOM=photomstack(f,fwhm,nom,2.5,3,6)
    photom.append(PHOTOM[0][0][1])
  return photom
  
"""
for i in range(6):
  f=freq[i]
  fwhm=psf(1.5)[i] #1.5= taille du pixel
  RES_PHOTOM.append(photomstack(f,fwhm)) 
  photom.append(RES_PHOTOM[i][0][0][1])
  noise.append(RES_PHOTOM[i][0][0][2])
"""
  
def th_KCMB(freq, A):
  x=freq/56.8
  return A*(x/(np.tanh(x/2))-4)

def fit_th_KCMB(donnees,freq): 
  data=donnees[:-1]
  freq=freq[:-1]
  popt,pcov=curve_fit(th_KCMB,freq,data,p0=(0.002)) #A initial 0.002
  f=np.arange(freq[0],freq[4])
  diff=data-th_KCMB(freq,*popt)
  plt.figure()
  plt.subplot(1,2,1)
  plt.scatter(freq, data, label='donnees', color='c')
  plt.plot(f, th_KCMB(f,*popt), label='fit KCMB, A=%5.3f'%tuple(popt))
  plt.plot(f,f*0,'k:')
  plt.legend()
  plt.title(r'Fit theorique $Flux=A\frac{x}{tanh(x/2)}-4$ avec x=f/56.8')
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique')
  plt.subplot(1,2,2)
  plt.scatter(freq,diff, color='g', marker='*')
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique')
  plt.title('Difference donnees-fit')
  plt.show()
  return popt,pcov

def th_MJY(freq, A):
  x=freq/56.8
  Flux_KCMB=x/(np.tanh(x/2))-4
  return A*Flux_KCMB*x**4*np.exp(x)/(np.exp(x)-1)**2

def fit_th_MJY(donnees,freq): 
  data=donnees[:-2]
  freq=freq[:-2]
  popt,pcov=curve_fit(th_MJY,freq,data,p0=(0.002)) #A initial 0.002
  f=np.arange(freq[0],freq[3])
  diff=data-th_MJY(freq,*popt)
  plt.figure()
  plt.subplot(1,2,1)
  plt.scatter(freq, data, label='donnees', color='c')
  plt.plot(f, th_MJY(f,*popt), label='fit MJY, A=%5.3f'%tuple(popt))
  plt.plot(f,f*0,'k:')
  plt.legend()
  plt.title(r'Fit theorique $Flux=AxF_{CMB}\frac{x^4e^x}{(e^x-1)^2}$ avec x=f/56.8')
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique')
  plt.subplot(1,2,2)
  plt.scatter(freq,diff, color='g', marker='*')
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique')
  plt.title('Difference donnees-fit')
  plt.show()
  return popt,pcov

def dust(freq,A,beta):
  freq=freq*10**9
  freq0=100e9
  #beta=1.8
  T=18
  k=1.38e-23
  h=6.63e-34
  return A*1e-6*(freq/freq0)**(beta+3)/(np.exp(h*freq/(k*T))-1) #1e-6 pour MJy

def dust_SZ_MJY(freq,Ad,beta,Asz):
  return dust(freq,Ad,beta)+th_MJY(freq, Asz)

def fit_dust_SZ_MJY(data,freq,r3): 
  popt,pcov=curve_fit(dust_SZ_MJY,freq,data,p0=(1,1.8,0.001)) #A initial 0.002
  f=np.arange(freq[0],freq[5])
  diff=data-dust_SZ_MJY(freq,*popt)
  plt.figure()
  plt.subplot(1,2,1)
  #plt.scatter(freq, data, label='donnees', color='c') #changer
  plt.errorbar(freq, data, statboot(freq,'MJY',r3)[1], label='donnees', color='c', fmt='o')
  plt.plot(f, dust_SZ_MJY(f,*popt), label='fit MJY, A_dust=%5.3f,beta=%5.3f,A_sz=%5.3f'%tuple(popt))
  plt.plot(f,f*0,'k:')
  plt.plot(f,dust(f,popt[0], popt[1]),'y--', label='Fit emission poussiere type corps gris')
  plt.plot(f,th_MJY(f,popt[2]),'g--', label='Fit emission SZ theorique')
  plt.legend()
  plt.title("Fit theorique emission SZ et emission des poussieres dans l'amas")
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique (MJY/sr)')
  plt.subplot(1,2,2)
  plt.errorbar(freq, diff, statboot(freq,'MJY',r3)[1], color='g', fmt='o')
  plt.plot(f,f*0,'k:')
  plt.ylim(-1,1)
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique (MJY/sr)')
  plt.title('Difference donnees-fit')
  plt.show()
  return popt,pcov

"""
plt.figure()
plt.title('Photometrie en fonction de la frequence')
plt.errorbar(freq,photom,noise, fmt='o', color='r', label='donnees')
plt.scatter(freq, Flux_KCMB, label='theorie KCMB') #*0.002
#plt.scatter(freq, Flux_JY, label='theorie mJY')
plt.axhline(y=0,xmin=0,xmax=900)
plt.legend()
plt.show()
"""

def plotphotom(f): #indice de la frequence
  plt.figure()
  source=RES_PHOTOM6[f][2]*RES_PHOTOM6[f][4]
  couronne=RES_PHOTOM6[f][3]*RES_PHOTOM6[f][4]
  data=RES_PHOTOM6[f][4]
  plt.subplot(1,3,1)
  plt.imshow(source.T,vmin=np.min(data),vmax=np.max(data),origin='lower')
  plt.subplot(1,3,2)
  plt.imshow(couronne.T,vmin=np.min(data),vmax=np.max(data),origin='lower')
  plt.subplot(1,3,3)
  plt.imshow(data.T,vmin=np.min(data),vmax=np.max(data),origin='lower')
  plt.show()



def photomboot(f,freq,nom, Nb,r3): #f la frequence et freq tableau des frequences
  stackboot=fits.open('../data/bootstrap_'+nom+'_'+str(f)+'_'+str(Nb)+'.fits')[0].data # changer
  photometrie=[]
  fwhm=psf(1.5)[np.where(freq==f)][0]
  for i in range(Nb):
    stack=stackboot[i,:,:]
    shape=np.shape(stack)
    peaks = np.zeros(shape, dtype=bool)
    peaks[64,64]=True
    photom = aperture_photom(stack, peaks, fwhm, 2.5, 3, r3)
    photometrie.append(photom[0][0][1])
  writefits(photometrie, 'Photometrie_bootstrap_'+nom+'_'+str(f)+'_'+str(Nb)+'_'+str(r3)+'.fits')
    

def statboot(freq,nom,r3): #freq tableau des frequences
  moyenne=[]
  ecartype=[]
  for i in freq:
    boot=fits.open('../data/Photometrie_bootstrap_'+nom+'_'+str(i)+'_1000_'+str(r3)+'.fits')[0].data
    moyenne.append(np.mean(boot))
    ecartype.append(np.std(boot))
  moyenne=np.array(moyenne)
  ecartype=np.array(ecartype)
  return moyenne, ecartype


#comparer nos moyennes a celles du bootstrap

def datavsboot(data,boot): #data:photometrie donnees boot:moyenne et erreur boot
  plt.figure()
  plt.errorbar(freq, boot[0], boot[1], label='bootstrap', color='c', fmt='o')
  plt.scatter(freq, data, color='g', marker='*', label='donnes')
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique (mJY/sr)')
  plt.legend()
  plt.title(r'Moyenne des flux photométriques du bootstrap (1000 stacks) VS flux photométriques pour nos données en fonction de la fréquence')
  plt.show()

#Etude ouvertures

def ouvertures(indf): #fonction permettant de determiner rayon couronne pour differentes ouvertures
  fwhm=psf(1.5)[indf]
  FR1=np.arange(1,3,0.1)
  FR2=[]
  FR3=[]
  frame_size = 8*int(fwhm)
  for i in FR1:
    r1=i*fwhm
    disk = round_mask((frame_size, frame_size), r1)
    r2=r1+fwhm/2
    disk_area=np.sum(disk)
    r3=np.sqrt((disk_area+np.pi*r2**2)/np.pi)
    FR2.append(r2/fwhm)
    FR3.append(r3/fwhm)
  return (FR1,FR2,FR3)

def photom_ouv(indf):
  freq=np.array([100,143,217,353,545,857])
  f=freq[indf]
  O=ouvertures(indf)
  FR1=O[0]
  FR2=O[1]
  FR3=O[2]
  Flux=[]
  stack=fits.open('../data/Stack_MJY_confirmed_amas_'+str(f)+'.fits')[0].data # changer
  shape=np.shape(stack)
  peaks = np.zeros(shape, dtype=bool)
  peaks[64,64]=True
  fwhm=psf(1.5)[indf]
  for i in range(len(FR1)):
    fr1=FR1[i]
    fr2=FR2[i]
    fr3=FR3[i]
    Flux.append(aperture_photom(stack, peaks, fwhm, fr1,fr2,fr3)[0][0][1])
  return Flux, FR1



"""
plt.figure()
plt.scatter(photom_ouv(0)[1],photom_ouv(0)[0])
plt.show()
"""
#plot des differents choix de photometrie


def fit_MJY_photom(data1,data2,data3,freq): 
  popt1,pcov1=curve_fit(dust_SZ_MJY,freq,data1,p0=(1,1.8,0.001))
  popt2,pcov2=curve_fit(dust_SZ_MJY,freq,data2,p0=(1,1.8,0.001))
  popt3,pcov3=curve_fit(dust_SZ_MJY,freq,data3,p0=(1,1.8,0.001))
  f=np.arange(freq[0],freq[5])
  plt.figure()
  #plt.scatter(freq, data, label='donnees', color='c') #changer
  plt.errorbar(freq, data1, statboot(freq,'MJY',4)[1], label='r1=2.5, r2=3, r3=4', color='c', fmt='o')
  plt.plot(f, dust_SZ_MJY(f,*popt1),'c--', label='fit mJY, A_dust=%5.3f,beta=%5.3f,A_sz=%5.3f'%tuple(popt1))
  plt.plot(f,f*0,'k:')
  plt.errorbar(freq, data2, statboot(freq,'MJY',5)[1], label='r1=2.5, r2=3, r3=5', color='g', fmt='o')
  plt.plot(f, dust_SZ_MJY(f,*popt2),'g--', label='fit mJY, A_dust=%5.3f,beta=%5.3f,A_sz=%5.3f'%tuple(popt2))
  plt.errorbar(freq, data3, statboot(freq,'MJY',6)[1], label='r1=2.5, r2=3, r3=6', color='y', fmt='o')
  plt.plot(f, dust_SZ_MJY(f,*popt3),'y--', label='fit mJY, A_dust=%5.3f,beta=%5.3f,A_sz=%5.3f'%tuple(popt3))
  plt.legend()
  plt.title("Fit theorique emission SZ et emission des poussieres dans l'amas pour differents choix d'ouvertures")
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique (mJY/sr)')
  plt.show()


def fit_MJY_cand_conf(confirme,candidat,freq):  #confirme photom6 cand=photometrie nom 
  popt1,pcov1=curve_fit(dust_SZ_MJY,freq,confirme,p0=(1,1.8,0.001)) #photom6
  #popt2,pcov2=curve_fit(dust_SZ_MJY,freq,candidat,p0=(1,1.8,0.001))
  f=np.arange(freq[0],freq[5])
  plt.figure()
  plt.errorbar(freq, confirme, statboot(freq,'MJY',6)[1], label='1094 Amas confirmes REDSHIFT>0', color='c', fmt='o')
  plt.plot(f, dust_SZ_MJY(f,*popt1),'c--', label='fit MJY, A_dust=%5.3f,beta=%5.3f,A_sz=%5.3f'%tuple(popt1))
  plt.plot(f,f*0,'k:')
  plt.errorbar(freq, candidat, statboot(freq,'MJY_candidats_amas',6)[1], label='1447 Candidats amas avec Q_NEURAL>0.7', color='g', fmt='o')
  #plt.plot(f, dust_SZ_MJY(f,*popt2),'g--', label='fit MJY, A_dust=%5.3f,beta=%5.3f,A_sz=%5.3f'%tuple(popt2))
  plt.legend()
  plt.title("Fit theorique emission SZ et emission des poussieres pour amas confirmes et candidats avec Q_NEURAL>0.7")
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique (MJY/sr)')
  plt.show()

def plot_fluxphotom(flux,freq):  #confirme photom6 cand=photometrie nom 
  f=np.arange(freq[0],freq[5])
  plt.figure()
  plt.scatter(freq, flux)
  plt.plot(f,f*0,'k:')
  plt.title("Flux photométrique en fonction de la fréquence pour 44 candidats (Hubert et Stephano)")
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique (MJY/sr)')
  plt.show()


def photometrie_bins(nom,freq): #nom MJY_bins_q
  PHOTOM_BINS=[]
  for i in range(6):
    f=freq[i]
    fwhm=psf(1.5)[i]
    stack_bins=fits.open('../data/Stack_MJY_bins_q_'+str(f)+'.fits')[0].data
    PHOTOM=[]
    for j in range(0,10):
      coupe=stack_bins[j,:,:]
      shape=np.shape(coupe)
      peaks = np.zeros(shape, dtype=bool)
      peaks[64,64]=True
      photom=aperture_photom(coupe,peaks,fwhm,2.5,3,6)
      PHOTOM.append(photom[0][0][1])
    PHOTOM_BINS.append(PHOTOM)
  PHOTOM_BINS=np.array(PHOTOM_BINS)
  return PHOTOM_BINS #pour toutes les frequences donne la photometrie de 10 stacks


def plot_q_bins(photom_bins,freq):  
  f=np.arange(freq[0],freq[5])
  plt.figure()
  plt.scatter(freq, photom_bins[:,0], label='0-165 Q_NEURAL', color='b')
  plt.scatter(freq, photom_bins[:,1], label='165-330 Q_NEURAL', color='g')
  plt.scatter(freq, photom_bins[:,2], label='330-495 Q_NEURAL', color='y')
  plt.scatter(freq, photom_bins[:,2], label='495-660 Q_NEURAL', color='r')
  plt.scatter(freq, photom_bins[:,3], label='660-825 Q_NEURAL', color='orange')
  plt.scatter(freq, photom_bins[:,4], label='825-990 Q_NEURAL', color='k')
  plt.scatter(freq, photom_bins[:,5], label='990-1155 Q_NEURAL', color='c')
  plt.scatter(freq, photom_bins[:,6], label='1155-1320 Q_NEURAL', color='xkcd:light pink')
  plt.scatter(freq, photom_bins[:,7], label='1320-1485 Q_NEURAL', color='m')
  plt.scatter(freq, photom_bins[:,8], label='1485-1650 Q_NEURAL', color='violet')
  plt.plot(f,f*0,'k:')
  plt.legend()
  plt.title("Flux photometrique pour differents bins Q_NEURAL")
  plt.xlabel('Frequence (GHz)')
  plt.ylabel('Flux photometrique (MJY/sr)')
  plt.show()
