#The main cleaning code. Requires certain files to sit
#in the same directory as it, and to specify certin initial
#parameters, as well as some working and output directories.

#In the same directory as this code, you need:
#-'triangle.py'
#-'binning.py'
#-'ccd4mask2.fit'
#-'iphas-images147_2_best.fits'
#-Your list of bright stars for the region from 'stars2.py'
#In the supplpied case, the list is called 'S147Stars0-14.fits'
#-'iphas-qc.fits'
#-In the same directory as this code there must contain
#the folder 'airglow_frames' containing the dark-time
#frames for all 3 filters.
#-The directory containing this code must also contain 
#a folder called 'confmaps'.
#-The root folder where the directory structure will be 
#built for all outputs must be set.



#moduls that need to be installed in Python:
#matplotlib
#math
#pylab
#numpy
#astropy
#collections
#urllib, urlparse
#aplpy
#PIL
#emcee
#scipy
#decimal
#pyfits


#to convert to laptop, only said filename3=filename2 and not to delete filename3
#removed funpack and rm .fz from image and confidence maps downloads
#added airglow subtraction from r-Ha image
#changed filenames for images and confidence maps
#changed root folder
import sys
import math
import matplotlib
import matplotlib.pyplot as plt
#plt.switch_backend('Qt4Agg') 
from pylab import *
import pylab
import numpy as np
from astropy.io import fits
import collections
import urllib, urlparse
import os
import aplpy
import binning as bin
from PIL import Image
#import suppress_stdout
import emcee
import triangle
from astropy import log
from scipy.interpolate.rbf import Rbf
from astropy.table import Table
from decimal import *
from astropy.wcs import WCS
import pyfits
from numpy import mean, sqrt, square, arange


#Set to begin and end at set image numbers in the list for testing,
#Must remove and replace commented out data when ready!!!

#####################INPUT-CONSTRAINTS########################

#Set name of run- will be top directory name containing outputs.
name = 'S147'

#Set confidence threshold for confidence map masking (0-100).
confidencethreshold=90

#Set size of bins, in pixels, for star removal.
stepi = 100 

#Set number of sigma to clip in input image.
sigLevel = 5 

#Set central galactic coordinates for area of sky to clean.
riascention = 180.2  	#l
declination = -1.6  	#b

#Set x and y sizes in degrees of area of sky to clean.
xsize = 1 		#x
ysize = 1 		#y

#Set root folder, in which to build directory structure and 
#for the code to work and store output files.
topdir = '/Users/charlesroe/Desktop/'

#Name (and location if not in same folder as this code) of 
#the fits table derivation of the original 'iphas-images.fits'.
fitstablebllocation = 'iphas-images.fits'#"iphas-images147_2_best.fits"

#Name (and location) of bright star table for the chosen region.
starTable_hdu = fits.open('IPHASStars0-16.fits')

##############################################################

root = topdir + name + '/'


#Opens tables and assigns the desired collumns to variables.
starTable = starTable_hdu[1].data
maxIt=1
iterations = np.arange(1,maxIt+1,1)
threshold = int(confidencethreshold)
iphas_qc = fits.open('iphas-qc.fits')[1].data
moonRuns = iphas_qc.field(76)
moon_phase = iphas_qc.field(71)
moon_altitude = iphas_qc.field(69)
moon_separation = iphas_qc.field(70)



#Opening arrays to store final data in output table.
one=[] 
two=[] 
three=[] 
four=[] 
five=[] 
six=[] 
seven=[] 
eight=[] 
nine=[] 
ten=[] 
eleven=[] 
twelve=[] 
thirteen=[] 
fourteen=[] 
fifteen=[] 
sixteen=[] 
seventeen=[] 
eighteen=[]
nineteen=[] 
twenty=[] 
twentyone=[] 
twentytwo=[] 
twentythree=[] 
moon_separations=[] 
exptimesr=[]
rtoHa_scalefactors=[] 
obsdatesHa=[] 
obsdatesr=[] 
juliandatesHa=[] 
juliandatesr=[] 
modjuliandatesHa=[]
modjuliandatesr=[] 
deltatees=[] 
rawHameds=[] 
rawrmeds=[] 
moon_phases=[] 
moon_altitudes=[]


#Setting coordinate system and names of filters.
coordsys, filters, filtr, filtHa = 'galactic', 'i', 'r', 'halpha'

#calculating the dimensions of the search area of sky.
armax = float(riascention)+float(xsize)/2
armin = float(riascention)-float(xsize)/2
admax = float(declination)+float(ysize)/2
admin = float(declination)-float(ysize)/2
size = str(xsize) + "x" + str(ysize)

#Opening iphas-images.fits
pointings_hdu = fits.open(fitstablebllocation)
pointings = pointings_hdu[1].data

#Building root directory structure.
os.system("mkdir " + root)
os.system("cd " + root + " && mkdir rawdir")
os.system("cd " + root + " && mkdir corrdir")
os.system("cd " + root + " && mkdir projdir")
os.system("cd " + root + " && mkdir diffdir")
os.system("cd " + root + " && mkdir final")
os.system("cd " + root + " && mkdir cleandir")
os.system("cd " + root + " && mkdir Ha_cleandir")
os.system("cd " + root + " && mkdir Colour_Plots")
os.system("cd " + root + " && mkdir Triangle_Plots")
os.system("cd " + root + " && mkdir models")
os.system("cd " + root + " && mkdir iterations")
os.system("cd " + root + " && mkdir MCMC_tables")
os.system("cd " + root + " && mkdir working_fol")
workfol = root+'working_fol/'
os.system("cd " + root + " && mkdir working_fol/Ha")
os.system("cd " + root + " && mkdir working_fol/r")
Hafol = root+'working_fol/Ha/'
rfol = root+'working_fol/r/'
os.system("cd " + root + " && mkdir working_fol/r_noHa")
rHafol = root+'working_fol/r_noHa/'
print 'root'
print 'Root built.'
print ' '

#Getting parameters from star table.
starRun = starTable.field(0)
starCCD = starTable.field(1)
starX = starTable.field(5)
starY = starTable.field(6)
starMag = starTable.field(4)
url = starTable.field(7)
starTable_hdu.close()

print 'Beginning search for extra-terrestrial life...'

#Get parameters from iphas-images.
l = pointings.field(24)
b = pointings.field(25)
in_dr2 = pointings.field(8)
best = pointings.field(9)
ra_min = pointings.field(26)
ra_max = pointings.field(28)
dec_min = pointings.field(27)
dec_max = pointings.field(29)
imagelist = pointings.field(2)
confmaplist = pointings.field(19)
filt = pointings.field(5)
ccd = pointings.field(1)
run = pointings.field(0)
field = pointings.field(7)
obstime = pointings.field(12)
pointings_hdu.close()

#Filter in_dr2 for bests only.
bests = best == True

print ' '
#in_dr2 = in_dr2 == 'true'
#in_dr2_2 = in_dr2[bests]

#Set up filter for only images in search area.
rcond1 = ra_min < armax
rcond2 = ra_max > armin
dcond1 = dec_min < admax
dcond2 = dec_max > admin
filts = filt == filters
filtsr = filt == filtr
filtsHa = filt == filtHa
#filts2 = filt == filters2
condsat = rcond1 & rcond2 & dcond1 & dcond2 & filts
condsatr = rcond1 & rcond2 & dcond1 & dcond2 & filtsr
condsatHa = rcond1 & rcond2 & dcond1 & dcond2 & filtsHa
#condsat_ha = rcond1 & rcond2 & dcond1 & dcond2 & filts2

#Filter all for bests only (in or out of dr2).
condsat2 = condsat[bests]
condsatr2 = condsatr[bests]
condsatHa2 = condsatHa[bests]
#condsat_ha2 = condsat_ha[bests]
imagelist2 = imagelist[bests]
confmaplist2 = confmaplist[bests]
l2 = l[bests]
b2 = b[bests]
run2 = run[bests]
ccd2 = ccd[bests]
field2=field[bests]
obstime2=obstime[bests]

##filter for only images in search area
#imagelistcond = imagelist2[condsat2]
#confmaplistcond = confmaplist2[condsat2]
#confmaplist_modified = []
#lcondsat = l2[condsat2]
#bcondsat = b2[condsat2]
#run = run2[condsat2]
#ccd = ccd2[condsat2]
#field=field2[condsat2]
#obstime=obstime2[condsat2]

#Apply area and filtermasks to arrays.
imagelistcondr = imagelist2[condsatr2]
confmaplistcondr = confmaplist2[condsatr2]
runr = run2[condsatr2]
ccdr = ccd2[condsatr2]
lcondsatr = l2[condsatr2]
bcondsatr = b2[condsatr2]
fieldr=field2[condsatr2]
obstimer=obstime2[condsatr2]
imagelistcondHa = imagelist2[condsatHa2]
confmaplistcondHa = confmaplist2[condsatHa2]
runHa = run2[condsatHa2]
ccdHa = ccd2[condsatHa2]
lcondsatHa = l2[condsatHa2]
bcondsatHa = b2[condsatHa2]
fieldHa=field2[condsatHa2]
obstimeHa=obstime2[condsatHa2]


if len(imagelistcondr) == len(imagelistcondHa):
	print 'No. of images: '+str(len(imagelistcondHa))
else:
	print 'No. of r images and Ha images do not match.'
	con = raw_input('Continue? (y/n): ')
	if con == 'n':
		quit()
print ' '

Begin = int(raw_input('Start at: '))
print ' '


#This loop downloads all the raw images in the search area into the specified directory.
sys.stdout.write("Downloading and Unpacking Raw Images")
if Begin == 0:
	notfixed = 0
	h = Begin
	while h<len(imagelistcondr):
	

		iir = runr[h]
		ccr = ccdr[h]
		
		iiHa = iir-1
		ccHa = ccr
		isit = iiHa in runHa
		sys.stdout.write(".")
		sys.stdout.flush()

		#print 'Downloading and Unpacking Raw Image....'+str(h+1)+'/'+str(len(imagelistcondr))
		if isit != True:
			print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
			print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
			print '!!!!!!!!!!!!!!!!!!!!!!!!!!!NOT RIGHT IMAGE!!!!!!!!!!!!!!!!!!!!!!!!!!'
			print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
			
			iiHa -= 1
			isit2 = iiHa in runHa
			if isit2 != True:
				print '!!!!!!!!!!!!!!!!!!!!!!!!!!STILL NOT CORRECT!!!!!!!!!!!!!!!!!!!!!!!!!'
				notfixed +=1
			elif isit2 == True:
				print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIXED WITH -2!!!!!!!!!!!!!!!!!!!!!!!!!!!'
				


		#imageName = 'r'+str(ii)+'-'+str(cc)
		imageNamer = 'r'+str(iir)+'-'+str(ccr)
		imageNameHa = 'r'+str(iiHa)+'-'+str(ccHa)
		#print imageNameHa, imageNamer
		
		#filename2 = str(imageName)+'_raw.fits.fz'
		filenamer2 = str(imageNamer)+'_raw.fits.fz'
		filenameHa2 = str(imageNameHa)+'_raw.fits.fz'

		filename3 = str(imageNamer)+'_raw.fits'
		filenameHa = str(imageNameHa)+'_raw.fits'

		#urllib.urlretrieve(str(imagelistcond[h]), workfol+filename2)
		url = 'http://www.iphas.org/data/images/'
		urlr = url+'r'+str(iir)[0:3]+'/'+imageNamer+'.fits.fz'
		urlHa = url+'r'+str(iiHa)[0:3]+'/'+imageNameHa+'.fits.fz'

		if str(iiHa) == '372903':
			print '!!!!!!!!!!!', urlr, urlHa


		urllib.urlretrieve(urlr, rfol+filename3)
		urllib.urlretrieve(urlHa, Hafol+filenameHa)

		#os.system('/home/croe/commands/cfitsio/funpack '+Hafol+str(filenameHa2))
		#os.system('/home/croe/commands/cfitsio/funpack '+rfol+str(filenamer2))
		#os.system('rm '+Hafol+str(filenameHa2))
		#os.system('rm '+rfol+str(filenamer2))
		h+=1
	print ' '
	print 'IMAGE DOWNLOADING COMPLETE.'
	print ' '
	print ' '
	print 'NO. OF IMAGES NOT FIXED WITH -2: '+str(notfixed)
print ' '



#The dark-time frames are opened.
airglow_r1 = pyfits.getdata('airglow_frames/fitted_r_1_airglow.fits')
airglow_r2 = pyfits.getdata('airglow_frames/fitted_r_2_airglow.fits')
airglow_r3 = pyfits.getdata('airglow_frames/fitted_r_3_airglow.fits')
airglow_r4 = pyfits.getdata('airglow_frames/fitted_r_4_airglow.fits')
airglow_Ha1 = pyfits.getdata('airglow_frames/fitted_Ha_1_airglow.fits')
airglow_Ha2 = pyfits.getdata('airglow_frames/fitted_Ha_2_airglow.fits')
airglow_Ha3 = pyfits.getdata('airglow_frames/fitted_Ha_3_airglow.fits')
airglow_Ha4 = pyfits.getdata('airglow_frames/fitted_Ha_4_airglow.fits')


print ' '
#The main cleaning loop begins.
ims = []
AS=[]
BS=[]
CS=[]
morethan1count=0
if Begin == 0:
	j = 0
else:
	j = Begin - 1
while j < len(imagelistcondr):#End:
	print str(j+1)+'/'+str(len(imagelistcondr))

#initial values set for iterative approach.
	totalModel = np.zeros((4096,2048))
	iterat = 0
	As=[]
	Bs=[]
	Cs=[]
	A=0
	B=0
	C=0


#Image name and CCD no. aquired.
	iir = runr[j]
	ccr = ccdr[j]		
	iiHa = iir-1
	ccHa = ccr
	isit = iiHa in runHa
	if isit != True:
		iiHa -= 1


#The information relating to the moon for this image is stored.
	moon_pha = float(moon_phase[moonRuns == int(iiHa)])
	moon_phases.append(moon_pha)
	moon_sep = float(moon_separation[moonRuns == int(iiHa)])
	moon_separations.append(moon_sep)
	moon_alt = float(moon_altitude[moonRuns == int(iiHa)])
	moon_altitudes.append(moon_alt)


#Observation length aquired.
	obsr_old = obstimer[j]




#Filenames for various steps are made.
	imageNamer = 'r'+str(iir)+'-'+str(ccr)
	imageNameHa = 'r'+str(iiHa)+'-'+str(ccHa)
	imageName=imageNamer
	filename = str(imageNamer)+'.fits'
	filenamer = str(imageNamer)+'.fits'
	filename2 = str(imageNamer)+'_raw.fits.fz'
	filename3 = str(imageNamer)+'_raw.fits'
	filename4 = str(imageNamer)+'_BS.fits'
	filenameHa = str(imageNameHa)+'_raw.fits'
	filenameHa2 = str(imageNameHa)+'_raw.fits.fz'
	filenameHa3 = str(imageNameHa)+'.fits'
	



	#filename3=filename2



#Raw images are opened and their observation lengths are grabbed from their headers.
	image_Ha = pyfits.getdata(Hafol+filenameHa)
	header_Ha = pyfits.getheader(Hafol+filenameHa, ext=1)
	obsHa = header_Ha['EXPTIME']
	image_r = pyfits.getdata(rfol+filename3)
	header_r = pyfits.getheader(rfol+filename3, ext=1)
	obsr = header_r['EXPTIME']
	#print header_Ha['OBJECT'], header_r['OBJECT']

#the image median is aquired.
	image_Ha_med = np.nanmedian(image_Ha)


#A load of header information is stored for the final output table.
	exptimesr.append(obsr)
	rawHameds.append(np.nanmedian(image_Ha))
	rawrmeds.append(np.nanmedian(image_r))
	obsdatesHa.append(header_Ha['DATE-OBS'])
	obsdatesr.append(header_r['DATE-OBS'])
	juliandatesHa.append(header_Ha['JD'])
	juliandatesr.append(header_r['JD'])
	modjuliandatesHa.append(header_Ha['MJD-OBS'])
	modjuliandatesr.append(header_r['MJD-OBS'])
	deltatee = float(header_r['JD'])-float(header_Ha['JD'])
	deltatees.append(deltatee)


#Scalefactors for the Ha and r data are calculated to bring them on to a comon exposure time.
	rHa_scalefactor = float(obsHa)/30. #(4)
	Ha_scalefactor = 30./float(obsHa) #(0.25)
	r_scalefactor = 30./float(obsr) #(1 or 3)
	#print Ha_scalefactor, r_scalefactor, obsHa, obsr, rHa_scalefactor


#Exposure time caled r-band image median.
	image_r_med = np.nanmedian(image_r*r_scalefactor)



#Finding median value of dark-time frames for dark-time scaling.
	if int(ccr) == 1:
		airglow_Ha_med = np.median(airglow_Ha1)
		airglow_r_med = np.median(airglow_r1)
	elif int(ccr) == 2:
		airglow_Ha_med = np.median(airglow_Ha2)
		airglow_r_med = np.median(airglow_r2)
	elif int(ccr) == 3:
		airglow_Ha_med = np.median(airglow_Ha3)
		airglow_r_med = np.median(airglow_r3)
	elif int(ccr) == 4:
		airglow_Ha_med = np.median(airglow_Ha4)
		airglow_r_med = np.median(airglow_r4)
	else:
		print '!!!!!!!![error 2]!!!!!!!'


#Limit of +4.0 above the dark-time median set. Any images with median below this 
#go no further in the cleaning process and have only their dark-time frames scaled a
#and subtracted.
	if image_Ha_med < (airglow_Ha_med + 4.):
		#subtract airglow frames (SCALE AIRGLOW FRAMES FIRST??)
		if int(ccr) == 1:

			airglow_scalefactor_Ha = image_Ha_med/airglow_Ha_med
			scaled_airglow_Ha = airglow_Ha1*airglow_scalefactor_Ha

		
			airglow_scalefactor_r = image_r_med/airglow_r_med
			scaled_airglow_r = airglow_r1*airglow_scalefactor_r


			image_r = (image_r*r_scalefactor) - scaled_airglow_r
			image_Ha = image_Ha - scaled_airglow_Ha
		
		elif int(ccr) == 2:

			airglow_scalefactor_Ha = image_Ha_med/airglow_Ha_med
			scaled_airglow_Ha = airglow_Ha2*airglow_scalefactor_Ha

		
			airglow_scalefactor_r = image_r_med/airglow_r_med
			scaled_airglow_r = airglow_r2*airglow_scalefactor_r


			image_r = (image_r*r_scalefactor) - scaled_airglow_r
			image_Ha = image_Ha - scaled_airglow_Ha


		elif int(ccr) == 3:

			airglow_scalefactor_Ha = image_Ha_med/airglow_Ha_med
			scaled_airglow_Ha = airglow_Ha3*airglow_scalefactor_Ha

		
			airglow_scalefactor_r = image_r_med/airglow_r_med
			scaled_airglow_r = airglow_r3*airglow_scalefactor_r


			image_r = (image_r*r_scalefactor) - scaled_airglow_r
			image_Ha = image_Ha - scaled_airglow_Ha


		elif int(ccr) == 4:

			airglow_scalefactor_Ha = image_Ha_med/airglow_Ha_med
			scaled_airglow_Ha = airglow_Ha4*airglow_scalefactor_Ha

		
			airglow_scalefactor_r = image_r_med/airglow_r_med
			scaled_airglow_r = airglow_r4*airglow_scalefactor_r


			image_r = (image_r*r_scalefactor) - scaled_airglow_r
			image_Ha = image_Ha - scaled_airglow_Ha


		else:
			print '!!!!!!!![error 2]!!!!!!!'

		print 'AIRGLOW SCALING: ', airglow_Ha_med, airglow_r_med, airglow_scalefactor_Ha, airglow_scalefactor_r
		#print np.nanmedian(image_Ha), np.nanmedian(image_r)	


#If an image median does not fall within the median of its corresponding dark-time 
#frame +4.0, then the dark-time frame is subtracted without any scaling and 
#the image continues to the next step.
	else:
		if int(ccr) == 1:
			image_r = (image_r*r_scalefactor) - airglow_r1
			image_Ha = image_Ha - airglow_Ha1

		elif int(ccr) == 2:
			image_r = (image_r*r_scalefactor) - airglow_r2
			image_Ha = image_Ha - airglow_Ha2

		elif int(ccr) == 3:
			image_r = (image_r*r_scalefactor) - airglow_r3
			image_Ha = image_Ha - airglow_Ha3

		elif int(ccr) == 4:
			image_r = (image_r*r_scalefactor) - airglow_r4
			image_Ha = image_Ha - airglow_Ha4

		else:
			print '!!!!!!!![error 2]!!!!!!!'





#Scale Ha to same exposure time as r (30 seconds).
	image_Ha = image_Ha*Ha_scalefactor
	

#r to Ha scalefactor after dark-time remoal calculated and stored.
	rtoHa_scalefactor = np.nanmedian(image_r)/np.nanmedian(image_Ha)
	if rtoHa_scalefactor < 12.00:
		rtoHa_scalefactor = 12.00
	rtoHa_scalefactors.append(rtoHa_scalefactor)


#Subtracting the Ha image from the r image to make an r-Ha image and saving.
	r_noHa = image_r-image_Ha
	pyfits.writeto(rHafol+filename3, r_noHa, header=header_r, clobber=True)



	


#Scaled r-Ha image opened.
	image_hdu = fits.open(rHafol+filename3)
	image =  np.array(image_hdu[0].data)
	image_header = image_hdu[0].header
	image_header1 = pyfits.getheader(rHafol+filename3, ext=0)

#URL set for confidence map download.
	url2 = 'http://www.iphas.org/data/images/confmaps/'

#Confidence map for r-band image downloaded.
	confnamer2 = confmaplistcondr[j]
	#print str(confname2)+ '!!!!!CONF!!!!!!'
	split = urlparse.urlsplit(str(confnamer2))
	aa = confnamer2.split('/')
	confnamer = aa[0]+'_'+aa[1]
	#change : in array name to _ 
	confnamer = list(confnamer)
	confnamer = np.array(confnamer)
	confnamer[confnamer==':'] = '_'
	confnamer = ''.join(confnamer)
	urllib.urlretrieve(str(url2)+str(confnamer2)+'.fz', 'confmaps/'+confnamer)#+'.fz')
	#os.system('/home/croe/commands/cfitsio/funpack confmaps/'+confnamer+'.fz')
	#os.system('rm confmaps/'+confnamer+'.fz')


#r-band confidence map opened.
	scidatar = pyfits.getdata('confmaps/'+str(confnamer), ext=int(ccr))


#Confidence map for Ha downloaded.
	confnameHa2 = confmaplistcondHa[j]
	#print str(confname2)+ '!!!!!CONF!!!!!!'
	split = urlparse.urlsplit(str(confnameHa2))
	aa = confnameHa2.split('/')
	confnameHa = aa[0]+'_'+aa[1]
	#change : in array name to _ 
	confnameHa = list(confnameHa)
	confnameHa = np.array(confnameHa)
	confnameHa[confnameHa==':'] = '_'
	confnameHa = ''.join(confnameHa)
	urllib.urlretrieve(str(url2)+str(confnameHa2)+'.fz', 'confmaps/'+confnameHa)#+'.fz')
	#os.system('/home/croe/commands/cfitsio/funpack confmaps/'+confnameHa+'.fz')
	#os.system('rm confmaps/'+confnameHa+'.fz')


#Ha confidence map opened.
	scidataHa = pyfits.getdata('confmaps/'+str(confnameHa), ext=int(ccHa))

	image_data = image

#Confidence maps masked to remove all pixels below confidence threshold.
	scidatar[scidatar<threshold] = 0
	scidatar[scidatar>=threshold] = 1
	scidatar = scidatar.astype(float)
	scidatar[scidatar == 0] = np.nan
	scidataHa[scidataHa<threshold] = 0
	scidataHa[scidataHa>=threshold] = 1
	scidataHa = scidataHa.astype(float)
	scidataHa[scidataHa == 0] = np.nan


	raw_header_Ha = pyfits.getheader(Hafol+filenameHa, ext=1)
	raw_image_Ha = pyfits.getdata(Hafol+filenameHa, ext=1)
	raw_image_r = pyfits.getdata(rfol+filename3, ext=1)
	raw_image_Ha_open = np.array([item for sublist in raw_image_Ha for item in sublist])





#Dark-time frames scaled and subtracted when image median is within the dark-time frame median + 4.0.
	if image_Ha_med < (airglow_Ha_med + 4.):
		#subtract airglow frames (SCALE AIRGLOW FRAMES FIRST??)
		if int(ccr) == 1:

			airglow_scalefactor_Ha = image_Ha_med/airglow_Ha_med
			scaled_airglow_Ha = airglow_Ha1*airglow_scalefactor_Ha

		
			airglow_scalefactor_r = image_r_med/airglow_r_med
			scaled_airglow_r = airglow_r1*airglow_scalefactor_r


			raw_image_r = (raw_image_r*r_scalefactor) - scaled_airglow_r
			raw_image_Ha = raw_image_Ha - scaled_airglow_Ha
		
		elif int(ccr) == 2:

			airglow_scalefactor_Ha = image_Ha_med/airglow_Ha_med
			scaled_airglow_Ha = airglow_Ha2*airglow_scalefactor_Ha

		
			airglow_scalefactor_r = image_r_med/airglow_r_med
			scaled_airglow_r = airglow_r2*airglow_scalefactor_r


			raw_image_r = (raw_image_r*r_scalefactor) - scaled_airglow_r
			raw_image_Ha = raw_image_Ha - scaled_airglow_Ha


		elif int(ccr) == 3:

			airglow_scalefactor_Ha = image_Ha_med/airglow_Ha_med
			scaled_airglow_Ha = airglow_Ha3*airglow_scalefactor_Ha

		
			airglow_scalefactor_r = image_r_med/airglow_r_med
			scaled_airglow_r = airglow_r3*airglow_scalefactor_r


			raw_image_r = (raw_image_r*r_scalefactor) - scaled_airglow_r
			raw_image_Ha = raw_image_Ha - scaled_airglow_Ha


		elif int(ccr) == 4:

			airglow_scalefactor_Ha = image_Ha_med/airglow_Ha_med
			scaled_airglow_Ha = airglow_Ha4*airglow_scalefactor_Ha

		
			airglow_scalefactor_r = image_r_med/airglow_r_med
			scaled_airglow_r = airglow_r4*airglow_scalefactor_r


			raw_image_r = (raw_image_r*r_scalefactor) - scaled_airglow_r
			raw_image_Ha = raw_image_Ha - scaled_airglow_Ha


		else:
			print '!!!!!!!![error 3]!!!!!!!'

		print 'AIRGLOW SCALING 2: ', airglow_Ha_med, airglow_r_med, airglow_scalefactor_Ha, airglow_scalefactor_r
		#print np.nanmedian(raw_image_Ha), np.nanmedian(raw_image_r)	


#If image median isn't within 4.0 of the dark-time frame median then the dark-time frame is 
#subtracted unscaled.
	else:
		if int(ccr) == 1:
			raw_image_r = (raw_image_r*r_scalefactor) - airglow_r1
			raw_image_Ha = raw_image_Ha - airglow_Ha1

		elif int(ccr) == 2:
			raw_image_r = (raw_image_r*r_scalefactor) - airglow_r2
			raw_image_Ha = raw_image_Ha - airglow_Ha2

		elif int(ccr) == 3:
			raw_image_r = (raw_image_r*r_scalefactor) - airglow_r3
			raw_image_Ha = raw_image_Ha - airglow_Ha3

		elif int(ccr) == 4:
			raw_image_r = (raw_image_r*r_scalefactor) - airglow_r4
			raw_image_Ha = raw_image_Ha - airglow_Ha4

		else:
			print '!!!!!!!![error 3]!!!!!!!'


	

	print ' '

#Confidence masking applied, including the custom mask for CCD 4.
	if int(ccr)==4:

		ccd4mask = pyfits.getdata('ccd4mask2.fit', ext=0)

		ccd4mask[ccd4mask == 0] = 1
		ccd4mask[ccd4mask > 1] = 0
		ccd4mask = ccd4mask.astype(float)
		ccd4mask[ccd4mask == 0] = np.nan

		raw_image_r_conf = raw_image_r * scidatar * ccd4mask
		raw_image_Ha_conf = raw_image_Ha * scidataHa * ccd4mask
		mask = image_data * scidatar * ccd4mask

	else:
		mask = image_data * scidatar
		raw_image_r_conf = raw_image_r * scidatar
		raw_image_Ha_conf = raw_image_Ha * scidataHa



#Ratio of good/bad pixels in confidence cleaned image computed and stored in final table.
	raw_image_Ha_conf_1D = np.array([item for sublist in raw_image_Ha_conf for item in sublist])
	jelly=0
	banana=0
	jam=0
	while jam<len(raw_image_Ha_conf_1D):
		peanut=raw_image_Ha_conf_1D[jam]
		if peanut == peanut:
			banana+=1
		jam+=1
	GoodPixRatio = (float(banana)/8388608.0)*100
	GoodPixRatio = round(GoodPixRatio, 2)
	fourteen.append(GoodPixRatio)


#Image medians recalculated after dark-time subtraction.
	afterairglow_med_Ha = np.nanmedian(raw_image_Ha_conf)
	afterairglow_med_r = np.nanmedian(raw_image_r_conf)

	leftoverfrommodelr = afterairglow_med_r/rtoHa_scalefactor
	leftoverfrommodelHa = afterairglow_med_Ha/rtoHa_scalefactor

	print 'LEFT OVER FROM MODEL (Ha, r) = ', leftoverfrommodelHa, leftoverfrommodelr


#If the image median is less than the dark-time median + 4.0 the process stops here, storing 
#the collected data in the final data table before moving on to the next image.
	if image_Ha_med < (airglow_Ha_med + 4.):
		print '!!!!!!!!!!!!!!BREAK!!!!!!!!!!!!!!'
		#print iiHa, ccHa, image_Ha_med, ' < ', airglow_Ha_med, ' +4.0'


		#raw_image_r_conf[raw_image_r_conf < 0.] = 0.
		#raw_image_Ha_conf[raw_image_Ha_conf < 0.] = 0.
		pyfits.writeto(root+'cleandir/'+str(imageName)+'.fits', raw_image_r_conf, header_r, clobber=True)
		pyfits.writeto(root+'Ha_cleandir/'+filenameHa3, raw_image_Ha_conf, header_Ha, clobber=True)



#Appending collected data for final data table.
		one.append(header_Ha['OBJECT'])
		two.append(iiHa)
		three.append(ccHa)
		four.append(header_Ha['WFFBAND'])
		five.append(obsHa)	
		six.append(np.nan)
		seven.append(np.nan)
		eight.append(np.nan)
		sixteen.append(np.nan)
		seventeen.append(np.nan)
		eighteen.append(np.nan)
		nineteen.append(np.nan)
		twenty.append(np.nan)
		twentyone.append(np.nan)
		nine.append(np.nan)
		ten.append(np.nan)
		eleven.append(np.nan)
		twelve.append(np.nan)
		thirteen.append(np.nan)
		twentytwo.append(np.nanmedian(raw_image_r_conf))
		twentythree.append(np.nanmedian(raw_image_Ha_conf))




		j+=1

		continue

	print iiHa, ccHa, image_Ha_med, ' > ', airglow_Ha_med, ' +4.0'
	print ' '







	pyfits.writeto(rfol+filename, raw_image_r_conf, image_header1, clobber=True)
	pyfits.writeto(Hafol+filenameHa3, raw_image_Ha_conf, raw_header_Ha, clobber=True)



#Bright star matching and removing.
	kk = 0
	shelob = 0
	starkk = []
	while kk < len(starRun):
		sRun = starRun[kk]
		sCCD = starCCD[kk]
	
		if str(sRun) == str(iir) and str(sCCD) == str(ccr):
			shelob+=1
			starkk.append(kk)
		kk+=1

	starkk = np.array(starkk)

#If there is a match:
	if shelob > 0:

		print '!!!!!!!!!!!!!!!!!'
		print 'NO. OF BRIGHT STARS TO BE REMOVED FROM IMAGE: '+ str(len(starkk))
		print '!!!!!!!!!!!!!!!!!'
		print 'Mag: R(pix):'
		starkkNo = 0
		while starkkNo<len(starkk):

			smag = starMag[starkk[starkkNo]]
			xStar = starX[starkk[starkkNo]]
			yStar = starY[starkk[starkkNo]]

			#finding WCS type and converting coords to pix coords
			w = WCS(image_header)

			px, py = w.wcs_world2pix(float(xStar), float(yStar), 1)


			#calculating mask axes
			centa, centb = py, px		#centre (y, x)(4096, 2048)
			radius = (8192/(float(smag)**2))+((1/float(smag))*1000)+100#(1/(float(smag)*1e-3))+(800/float(smag))
			sr = int(radius)				#radius
			print smag, sr
			if sr > 2000:
				sr=2000
			#creating circular mask
			sy, sx = np.ogrid[-centa:4096-centa, -centb:2048-centb]
			smask = sx*sx + sy*sy <= sr*sr

			#applying mask
			mask[smask] = np.nan
			starkkNo+=1

#Saving bright-star-masked image
	pyfits.writeto(root+'iterations/'+filename4, mask, image_header1, clobber=True)
	#os.system('rm '+str(filename3))


	print ' '
	
#Begin iteraive cleaning approach.
	while iterat<len(iterations):
		iteration = int(iterations[iterat])
		print '!!!!!!!!!!!!!!!!'
		print str(imageName), 'Iteration: '+str(iteration)
		print '!!!!!!!!!!!!!!!!'
		print ' '
#Image names depending on whether it's the first iteration or not.
		if iteration > 1:
			filename = str(imageName)+'_'+str(iteration-1)+'.fits'
		else:
			filename = filename4


#Image opened.
		image = pyfits.getdata(root+'iterations/'+filename4, ext=0)
		image_header = pyfits.getheader(root+'iterations/'+filename4, ext=0)

#Image information saved for final data table.
		ims.append(imageName)
		one.append(header_Ha['OBJECT'])
		two.append(iiHa)
		three.append(ccHa)
		four.append(header_Ha['WFFBAND'])
		five.append(obsHa)	


#Version saved for colour plot later. This image is deleted before the end.
		pyfits.writeto(workfol+str(imageName)+'_imagez.fits', mask, clobber=True)

#Size of bins set.
		stepi=int(stepi)
		y_range = np.arange(0,4096,1)
		x_range = np.arange(0,2048,1)

		########areas3
		
		image2 = np.array([item for sublist in image for item in sublist])


		
#Parameters for binning the image to remove stars are set. These include the pixel coordinates for
#each bin and length of each axis of the bined image.
		y_range3 = np.arange(100,4100,stepi)
		x_range3 = np.arange(100,2100,stepi)
		xLen=int(len(x_range3)-1)
		yLen=int(len(y_range3)-1)


#The binning loop, where the confidence , dark-time and bright star cleaned image is split into
#bins and the medioan value of each bin is taken to ensure the model is fitting only the background
#of the image.
		xCoords=[]
		yCoords=[]
		medVals=[]
		trial=0
		while trial<(len(y_range3)-1):

			y_min=y_range3[trial]
			y_max=y_range3[trial+1]
			yCoord=(y_min+y_max)/2
			jeer=0
			while jeer<(len(x_range3)-1):
				binVals=[]
				y_min2=y_min
				x_min=x_range3[jeer]
				x_max=x_range3[jeer+1]
				xCoord=(x_min+x_max)/2
				while y_min2<y_max:
					arr = image[y_min]
					xar = arr[x_min:x_max]
					binVals.append(xar)

					y_min2+=1

				binVals=np.array(binVals)
				binVals = np.array([item for sublist in binVals for item in sublist])
				binVals_noNaNs = binVals[binVals==binVals]
				
				if float(len(binVals_noNaNs))/float(len(binVals))<0.7:
					medVal = np.nan
				else:
					medVal = np.nanmedian(binVals)


				
				medVals.append(medVal)
				xCoords.append(xCoord)
				yCoords.append(yCoord)
				jeer+=1

			trial+=1


#Final arrays containing the median value of each bin and the coordinates of each bin.
		medVals=np.array(medVals)
		xCoords=np.array(xCoords)
		yCoords=np.array(yCoords)


#rms of binned image before model
		bef = []
		bbb=0
		while bbb<len(medVals):
			bef1 = medVals[bbb]
			bef1 = bef1 - np.nanmedian(medVals)
			if bef1 == bef1:
				bef1=bef1**2
				bef.append(bef1)
			bbb+=1
		bef=np.array(bef)
		meanbef=np.mean(bef)
		rmsBefore = math.sqrt(float(meanbef))
		nine.append(rmsBefore)


#standard deviation taken of the binned image and sigma clipping takes place.
		med_stddev = np.nanstd(medVals)
		med_med = np.nanmedian(medVals)
		print str(med_med), str(med_stddev), 'med and stddv, and sigma limits below:'
		lower_limit = (med_med-(float(sigLevel)*med_stddev))
		upper_limit = (med_med+(float(sigLevel)*med_stddev))
		print upper_limit, lower_limit
		medVals[medVals!=medVals] = 0.
		medVals[medVals<lower_limit] = np.nan
		medVals[medVals>upper_limit] = np.nan
		medVals[medVals==0.] = np.nan
		print ' '

#The binned image is saved.
		pyfits.writeto(workfol+imageName+'_binnedforMCMC.fits', medVals.reshape(yLen,xLen), clobber=True)

#The min and max values for colour plots calculated.
		itemsccc2=[]
		ccc2=0
		while ccc2<len(medVals):
			itemccc2=medVals[ccc2]
			if itemccc2==itemccc2:
				itemsccc2.append(itemccc2)
			ccc2+=1
		itemsccc2=np.array(itemsccc2)
		stddev2=np.std(itemsccc2)
		med2 = np.median(itemsccc2)
		itemsccc2[itemsccc2>(med2+2*stddev2)] = med2
		itemsccc2[itemsccc2<(med2-2*stddev2)] = med2
		minVal=min(itemsccc2)
		maxVal=max(itemsccc2)


#The binned values are saved to a fits table for use by the emcee model fitter. (median values for each bin
	#wth their x and y pixel coordinates)
		col1 = fits.Column(name='x', format='E', array=xCoords)
		col2 = fits.Column(name='y', format='E', array=yCoords)
		col3 = fits.Column(name='z', format='E', array=medVals)
		cols = fits.ColDefs([col1, col2, col3])
		tbhdu = fits.BinTableHDU.from_columns(cols)
		tbhdu.writeto(workfol+str(imageName)+'_table_'+str(iteration)+'.fits', clobber=True)
		#########

#The emcee sampler that fits the bst fit model to the binned data.
		def lnprob(theta, xvals, yvals, data, errors):
		    lp = lnprior(theta)
		    if not np.isfinite(lp):
		        return -np.inf
		    return lnlike(theta, xvals, yvals, data, errors)

		def writeout(sampler):
			flatchains = mcmc_results(sampler,param_keys).flatchain
			outdata = Table(flatchains.values(), names=flatchains.keys())
			table_filename=(root+'MCMC_tables/'+str(imageName)+'_MCMCresults_it'+str(iteration)+'.fits')
			Table.write(outdata,table_filename,overwrite=True)

		def mcmctriangle(ID):
			table_filename=(root+'MCMC_tables/'+str(imageName)+'_MCMCresults_it'+str(iteration)+'.fits')
			data = Table.read(table_filename)
			data_t = np.array([data[key] for key in param_keys]).transpose()
			triangle.corner(data_t,labels=param_labels,quantiles=[0.1587,0.5000,0.8413])
			plt.savefig(root+'Triangle_Plots/'+str(imageName)+'_triangle.pdf', clobber=True)
			plt.close()
			#plt.show()

		class mcmc_results:
			def __init__(self, sampler, param_keys):
				self.param_keys = param_keys
				self.sampler = sampler
				self.flatchain = self.get_flatchains()

			def get_flatchains(self):
				results = [self.sampler.flatchain[:,i] for i in range(len(param_keys))]
				return dict(zip(self.param_keys, results))

			def chain(self,walker):
				chains = [self.sampler.chain[walker][:,i] for i in range(len(param_keys))]
				return dict(zip(self.param_keys, chains))

		def lnprior(theta):
		    a, b, c = theta
		    if (-10000. < c < 10000. and -1. < b < 1. and -1. < a < 1.):
		        return 0.0
		    return -np.inf

		def lnlike(theta, xvals, yvals, data, errors):
		    a, b, c = theta
		    model = a*xvals + b*yvals + c
		    diff = data-model
		    invsigma2 = 1./errors**2
		    return -0.5*(np.sum(diff**2*invsigma2 - np.log(invsigma2)))


	#Set parameters for emcee and set starting positio
		ndim, nwalkers, nruns, burnruns = 3, 100, 1000, 200
		#only need to change nwalkers, nruns, burnruns

		param_keys = ['a','b','c']
		param_labels=[r'a',r'b',r'c']

	#Read in data: a row for each source counts value, a column for each background value
		datatable_hdu = fits.open(workfol+str(imageName)+'_table_'+str(iteration)+'.fits')
		datatable = datatable_hdu[1].data

#Input the data poitns, as well as their x and y coordinates from the data table.
		xvals = datatable.field(0)
		yvals = datatable.field(1)
		data = medVals #datatable.field(2)#medians3#
	
#Mask any NaN values as 0.
		medVals[medVals!=medVals] = 0.0

#Errors are colculated for the data and stored. Any masked NaN values from the previous step 
#are given huge errors so the sampler does not include them when searching for the best fit.
		LowerPercentile = np.percentile(medVals, 16)
		MidPercentile = np.percentile(medVals, 50)
		UpperPercentile = np.percentile(medVals, 84)
		eleven.append(LowerPercentile)
		twelve.append(MidPercentile)
		thirteen.append(UpperPercentile)
		error1 = (UpperPercentile-LowerPercentile)/2
		error2 = 9e99
		err = []
		char = 0
		while char<len(medVals):
			point = medVals[char]
			if point == 0.0:
				err.append(error2)
			elif point != 0.0:
				err.append(error1)
			char+=1
		errors = np.array(err)#(UpperPercentile-LowerPercentile)/2   #data*0. + 1.

#Create random starting values
		pos = [np.array([0.1,0.1,0.1]) + np.random.rand(ndim)/10 for i in range(nwalkers)]

#Set up the emcee sampler
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(xvals, yvals, data, errors))
		pos, prob, state = sampler.run_mcmc(pos, burnruns)
		sampler.reset()
		sampler.run_mcmc(pos, nruns)

		writeout(sampler)
		mcmctriangle(0)

#The MCMC output table is opened to get the best fit values.
		resultstable_hdu = fits.open(root+'MCMC_tables/'+str(imageName)+'_MCMCresults_it'+str(iteration)+'.fits')
		resultstable = resultstable_hdu[1].data

#The best fit values for A, B and C are taken from the tables, their upper and lower percentiles are found,
# and all results are saved for the final data table.
		Ayes=resultstable.field(0)
		Bees=resultstable.field(2)
		Cees=resultstable.field(1)
		LowerA = np.percentile(Ayes, 16)
		UpperA = np.percentile(Ayes, 84)
		LowerB = np.percentile(Bees, 16)
		UpperB = np.percentile(Bees, 84)
		LowerC = np.percentile(Cees, 16)
		UpperC = np.percentile(Cees, 84)
		A = np.median(Ayes)
		B = np.median(Bees)
		C = np.median(Cees)
		six.append(A)
		seven.append(B)
		eight.append(C)
		sixteen.append(LowerA)
		seventeen.append(UpperA)
		eighteen.append(LowerB)
		nineteen.append(UpperB)
		twenty.append(LowerC)
		twentyone.append(UpperC)
		print ' '
		print 'A -> '+str(A)+' B -> '+str(B)+' C -> '+str(C)
		print ' '


#Those NaN values masked to 0.0 for the emcee sampler are returned to nan values.
		medVals[medVals==0.0] = np.nan




# The model is computed on a 1x1 scale using A, B and C.
		zVals6=[]
		par6=0
		meansaft=0
		while par6<(len(y_range-1)):
			y6= y_range[par6]
			per6=0
			while per6<(len(x_range-1)):
				x6 = x_range[per6]
				Z6=0
				Z6 = (A*x6)+(B*y6)+C
				zVals6.append(Z6)
				#sqaft = Z6**2
				#meansaft += sqaft
		
				per6+=1	
			par6+=1
		zVals6 = np.array(zVals6)
		#print zVals6
		MedOfMeds = np.median(medVals)
		model4 = zVals6.reshape((4096,2048))


#The model is then added to the total iterative model.
		totalModel += model4
		#print totalModel

#The total model is then subtracted from the image data.
		##compared 4 is now the residuals
		compared4 = image-totalModel
		##compared5 is the residuals plus the median from the raw data
		##compared5 = compared4+MedOfMeds

		##COMPARED4 SHOULD BE FINAL IMAGE ARRAY!!!!



#The model is computed on a 100x100 scale.
		y_range4 = np.arange(100,4100,stepi)
		x_range4 = np.arange(100,2100,stepi)
		zVals=[]
		par=0
		while par<(len(y_range4)-1):
			y = y_range4[par]
			per=0
			while per<(len(x_range4)-1):
				x = x_range4[per]
				Z = A*x+B*y+C
				zVals.append(Z)
				per+=1	
			par+=1
		zVals = np.array(zVals)

	

#A version of the model is saved as a fits image for use in the colour plots. this image is deleted before the end.
		model = zVals.reshape((yLen,xLen))
		pyfits.writeto(workfol+str(imageName)+'_modelz.fits', model, clobber=True)

#A version of the binned data is saved as a fits image for use in the colour plots. this image is deleted before the end.
		medians4 = medVals.reshape((yLen,xLen))
		pyfits.writeto(workfol+str(imageName)+'_binnedz.fits', medians4, clobber=True)

#A version of the residuals is saved as a fits image for use in the colour plots. this image is deleted before the end.
		compared2 = medVals-zVals
		compared = compared2#+MedOfMeds
		compared = np.array(compared)
		compare = compared.reshape((yLen,xLen))
		pyfits.writeto(workfol+str(imageName)+'_fixedz.fits', compare, clobber=True)

		filenamer = str(imageNamer)+'.fits'
#WRITING FINAL CLEANED IMAGE FILE TO DIRECTORY (unbinned)
		if iteration == maxIt:
			image_to_clean = pyfits.getdata(rfol+filenamer)
			Ha_image_to_clean = pyfits.getdata(Hafol+filenameHa3)
			Ha_image_header = pyfits.getheader(Hafol+filenameHa3, ext=0)
			#print filenamer, filenameHa3

			#scaled_model = totalModel*(12.84/11.84)
			scaled_model = totalModel 
			#Ha_model = (scaled_model/12.84)*rHa_scalefactor
			Ha_model = (scaled_model/rtoHa_scalefactor)*rHa_scalefactor
			

			cleaned_image = image_to_clean - scaled_model
			cleaned_image = cleaned_image - leftoverfrommodelr
			twentytwo.append(np.nanmedian(cleaned_image))
			Ha_cleaned_image = Ha_image_to_clean - Ha_model
			Ha_cleaned_image = Ha_cleaned_image - leftoverfrommodelHa
			twentythree.append(np.nanmedian(Ha_cleaned_image))
			pyfits.writeto(root+'cleandir/'+str(imageName)+'.fits', cleaned_image, image_header, clobber=True)
			pyfits.writeto(root+'Ha_cleandir/'+filenameHa3, Ha_cleaned_image, Ha_image_header, clobber=True)

		else:
			pyfits.writeto(root+'iterations/'+str(imageName)+'_'+str(iteration)+'.fits', compared4, image_header, clobber=True)
		#os.system('rm '+str(imageName)+'.fits.fz')
		#os.system('rm '+str(imageName)+'_results.fits')
		#os.system('rm '+str(imageName)+'_table.fits')



#residual image binned to calculate rms after cleaning.
		medVals3=[]
		trial=0
		while trial<(len(y_range3)-1):
			y_min=y_range3[trial]
			y_max=y_range3[trial+1]
			jeer=0
			while jeer<(len(x_range3)-1):
				binVals=[]
				y_min2=y_min
				x_min=x_range3[jeer]
				x_max=x_range3[jeer+1]
				while y_min2<y_max:
					arr = Ha_cleaned_image[y_min]
					xar = arr[x_min:x_max]
					binVals.append(xar)
					y_min2+=1
				binVals=np.array(binVals)
				binVals = np.array([item for sublist in binVals for item in sublist])
				binVals_noNaNs = binVals[binVals==binVals]
				if float(len(binVals_noNaNs))/float(len(binVals))<0.7:
					medVal = np.nan
				else:
					medVal = np.nanmedian(binVals)
				medVals3.append(medVal)
				jeer+=1
			trial+=1
		medVals3=np.array(medVals3)






#rms of binned residuals scaled
		aft = []
		aaa=0
		while aaa<len(medVals3):
			aft1 = medVals3[aaa]
			if aft1 == aft1:
				aft1=aft1**2
				aft.append(aft1)
			aaa+=1
		aft=np.array(aft)

		meanaft=np.mean(aft)
		rmsAfter=math.sqrt(float(meanaft))

		ten.append(rmsAfter)

##############################################################################
##############################################################################
#This whole next secion is dedicated to plotting the colour plots.
		colourMin=minVal
		colourMax=maxVal
		print 'vmin='+str(minVal),' vmax='+str(maxVal)
		print ' '
#The first panel, the original image.
		fig = plt.figure()
		fig.patch.set_facecolor('white')
		f1 = aplpy.FITSFigure(workfol+str(imageName)+'_imagez.fits', hdu=0, figure=fig, subplot=(1,5,1))
		f1.show_colorscale(vmin=colourMin, vmax=colourMax)
		f1.axis_labels.set_xtext('Pixels')
		f1.axis_labels.set_ytext('Pixels')
		f1.tick_labels.set_font(size='small', weight='medium', \
    	                     stretch='normal', family='sans-serif', \
    	                     style='normal', variant='normal')
		f1.show_grid()
		getcontext().prec=3
		f1.set_title('Original Image')



#The second panel, the binned image.
		f2 = aplpy.FITSFigure(workfol+str(imageName)+'_binnedz.fits', hdu=0, figure=fig, subplot=(1,5,2))
		f2.show_colorscale(vmin=colourMin, vmax=colourMax)
		f2.axis_labels.set_xtext('x100 Pixels')
		f2.hide_yaxis_label()
		f2.tick_labels.hide_y()
		f2.tick_labels.set_font(size='small', weight='medium', \
    	                     stretch='normal', family='sans-serif', \
    	                     style='normal', variant='normal')
		f2.show_grid()
		f2.set_title('Binned')


#the third panel, the model.
		f3 = aplpy.FITSFigure(workfol+str(imageName)+'_modelz.fits', hdu=0, figure=fig, subplot=(1,5,3))
		f3.show_colorscale(vmin=colourMin, vmax=colourMax)
		f3.axis_labels.set_xtext('x100 Pixels')
		f3.hide_yaxis_label()
		f3.tick_labels.hide_y()
		f3.tick_labels.set_font(size='small', weight='medium', \
    	                     stretch='normal', family='sans-serif', \
    	                     style='normal', variant='normal')
		f3.show_grid()
		f3.set_title('Model')



#The colourbar for the first 3 panels is set.
		cbar_ax1 = fig.add_axes([0.6, 0.23, 0.02, 0.485])
		cb1 = plt.colorbar(f2.image, cax = cbar_ax1)
		cb1.set_label('Pixel Count')



#The fourth panel, the residuals.
		f4 = aplpy.FITSFigure(workfol+str(imageName)+'_fixedz.fits', hdu=0, figure=fig, subplot=(1,5,5))
		f4.show_colorscale()
		f4.axis_labels.set_xtext('x100 Pixels')
		f4.hide_yaxis_label()
		f4.tick_labels.hide_y()
		f4.tick_labels.set_font(size='small', weight='medium', \
    	                     stretch='normal', family='sans-serif', \
    	                     style='normal', variant='normal')
		f4.show_grid()
		f4.set_title('Residuals')


#the colourbar for the fourth panel is set.
		cbar_ax2 = fig.add_axes([0.92, 0.23, 0.02, 0.485])
		cb2 = plt.colorbar(f4.image, cax = cbar_ax2)
		cb2.set_label('Pixel Count')


#the plot is saved.
		plt.savefig(root+'Colour_Plots/'+str(imageName)+'_plot'+str(iteration)+'.png', bbox_inches='tight', clobber=True)

#The images that make up each panel are deleted.
		os.system('rm '+workfol+str(imageName)+'_imagez.fits')
		os.system('rm '+workfol+str(imageName)+'_binnedz.fits')
		os.system('rm '+workfol+str(imageName)+'_modelz.fits')
		os.system('rm '+workfol+str(imageName)+'_fixedz.fits')
		plt.close()
		
##################################################################################
##################################################################################
	

	

#Various closures.
		plt.close()
		As.append(A)
		Bs.append(B)
		Cs.append(C)
		##delete previous iterations
		#if iteration > 1:
		#	os.system('''rm str(imageName)+'_'+str(iteration-1)+'.fits''') 
		image_hdu.close()
		datatable_hdu.close()
		resultstable_hdu.close()


		iterat+=1

#The total model is saved.
	pyfits.writeto(root+'models/'+str(imageName)+'_total_r_Model.fits', scaled_model, clobber=True)


#the evolution of A, B and C over iterations is plotted.
	As = np.array(As)
	Bs = np.array(Bs)
	Cs = np.array(Cs)
	As = np.absolute(As)
	Bs = np.absolute(Bs)
	Cs = np.absolute(Cs)
	As = (As/As[0])*100
	Bs = (Bs/Bs[0])*100
	Cs = (Cs/Cs[0])*100
	plt.plot(As, 'r-', label='A')
	plt.plot(Bs, 'g-', label='B')
	plt.plot(Cs, 'b-', label='C')
	plt.title(str(imageName) + ' Iterative Cleaning.')
	plt.xlabel('Iteration')
	plt.ylabel('Percentage of initial fit')
	plt.ylim((0,100))
	plt.legend(prop={'size':11})
	plt.savefig(root+'Colour_Plots/'+str(imageName)+'_Alls.png')
	plt.close()




#The final large fits data table containing all the calculated parameters for each image is stored after each image.
	col1 = fits.Column(name='field', format='20A', array=np.array(one))
	col2 = fits.Column(name='run', format='20A', array=np.array(two))
	col3 = fits.Column(name='ccd', format='20A', array=np.array(three))
	col4 = fits.Column(name='filter', format='20A', array=np.array(four))
	col5 = fits.Column(name='Ha_exp_time', format='E', array=np.array(five))
	col6 = fits.Column(name='A', format='E', array=np.array(six))
	col7 = fits.Column(name='B', format='E', array=np.array(seven))
	col8 = fits.Column(name='C', format='E', array=np.array(eight))
	col9 = fits.Column(name='rms_bef', format='E', array=np.array(nine))
	col10 = fits.Column(name='rms_aft', format='E', array=np.array(ten))
	col11 = fits.Column(name='Binned_low', format='E', array=np.array(eleven))
	col12 = fits.Column(name='Binned_med', format='E', array=np.array(twelve))
	col13 = fits.Column(name='Binned_upp', format='E', array=np.array(thirteen))
	col14 = fits.Column(name='good_pix_perc', format='20A', array=np.array(fourteen))
	col15 = fits.Column(name='comments', format='20A', array=np.array(fifteen))
	col16 = fits.Column(name='A_low', format='E', array=np.array(sixteen))
	col17 = fits.Column(name='A_upp', format='E', array=np.array(seventeen))
	col18 = fits.Column(name='B_low', format='E', array=np.array(eighteen))
	col19 = fits.Column(name='B_upp', format='E', array=np.array(nineteen))
	col20 = fits.Column(name='C_low', format='E', array=np.array(twenty))
	col21 = fits.Column(name='C_upp', format='E', array=np.array(twentyone))
	col22 = fits.Column(name='r_med_resid', format='E', array=np.array(twentytwo))
	col23 = fits.Column(name='Ha_med_resid', format='E', array=np.array(twentythree))
	col24 = fits.Column(name='r to Ha', format='E', array=np.array(rtoHa_scalefactors))
	col25 = fits.Column(name='Ha obs date', format='30A', array=np.array(obsdatesHa))
	col26 = fits.Column(name='r obs date', format='30A', array=np.array(obsdatesr))
	col27 = fits.Column(name='JD_Ha', format='20E', array=np.array(juliandatesHa))
	col28 = fits.Column(name='JD_r', format='20E', array=np.array(juliandatesr))
	col29 = fits.Column(name='modified_JD_Ha', format='20E', array=np.array(modjuliandatesHa))
	col30 = fits.Column(name='modified_JD_r', format='20E', array=np.array(modjuliandatesr))
	col31 = fits.Column(name='delta_t_(r-Ha)', format='10E', array=np.array(deltatees))
	col32 = fits.Column(name='raw_Ha_med', format='E', array=np.array(rawHameds))
	col33 = fits.Column(name='raw_r_med', format='E', array=np.array(rawrmeds))
	col34 = fits.Column(name='moon_altitude', format='E', array=np.array(moon_altitudes))
	col35 = fits.Column(name='moon_separation', format='E', array=np.array(moon_separations))
	col36 = fits.Column(name='moon_phase', format='E', array=np.array(moon_phases))
	col37 = fits.Column(name='r_exp_time', format='E', array=np.array(exptimesr))
	cols = fits.ColDefs([col1, col2, col3, col4, col5, col37, col32, col33, col24, col22, col23, col6, col7, col8, col9, col10, col11, col12, col13, col14, col16, col17, col18, col19, col20, col21, col25, col26, col27, col28, col29, col30, col31, col34, col35, col36, col15])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(root + 'iphas-cleaning_'+str(Begin)+'.fits', clobber=True)

#The end of the cleaning loop.
	j+=1



#The final version of the data table, after all images have been processed is stored.
print 'Final table store..'
col1 = fits.Column(name='field', format='20A', array=np.array(one))
col2 = fits.Column(name='run', format='20A', array=np.array(two))
col3 = fits.Column(name='ccd', format='20A', array=np.array(three))
col4 = fits.Column(name='filter', format='20A', array=np.array(four))
col5 = fits.Column(name='Ha_exp_time', format='E', array=np.array(five))
col6 = fits.Column(name='A', format='E', array=np.array(six))
col7 = fits.Column(name='B', format='E', array=np.array(seven))
col8 = fits.Column(name='C', format='E', array=np.array(eight))
col9 = fits.Column(name='rms_bef', format='E', array=np.array(nine))
col10 = fits.Column(name='rms_aft', format='E', array=np.array(ten))
col11 = fits.Column(name='Binned_low', format='E', array=np.array(eleven))
col12 = fits.Column(name='Binned_med', format='E', array=np.array(twelve))
col13 = fits.Column(name='Binned_upp', format='E', array=np.array(thirteen))
col14 = fits.Column(name='good_pix_perc', format='20A', array=np.array(fourteen))
col15 = fits.Column(name='comments', format='20A', array=np.array(fifteen))
col16 = fits.Column(name='A_low', format='E', array=np.array(sixteen))
col17 = fits.Column(name='A_upp', format='E', array=np.array(seventeen))
col18 = fits.Column(name='B_low', format='E', array=np.array(eighteen))
col19 = fits.Column(name='B_upp', format='E', array=np.array(nineteen))
col20 = fits.Column(name='C_low', format='E', array=np.array(twenty))
col21 = fits.Column(name='C_upp', format='E', array=np.array(twentyone))
col22 = fits.Column(name='r_med_resid', format='E', array=np.array(twentytwo))
col23 = fits.Column(name='Ha_med_resid', format='E', array=np.array(twentythree))
col24 = fits.Column(name='r to Ha', format='E', array=np.array(rtoHa_scalefactors))
col25 = fits.Column(name='Ha obs date', format='30A', array=np.array(obsdatesHa))
col26 = fits.Column(name='r obs date', format='30A', array=np.array(obsdatesr))
col27 = fits.Column(name='JD_Ha', format='20E', array=np.array(juliandatesHa))
col28 = fits.Column(name='JD_r', format='20E', array=np.array(juliandatesr))
col29 = fits.Column(name='modified_JD_Ha', format='20E', array=np.array(modjuliandatesHa))
col30 = fits.Column(name='modified_JD_r', format='20E', array=np.array(modjuliandatesr))
col31 = fits.Column(name='delta_t_(r-Ha)', format='10E', array=np.array(deltatees))
col32 = fits.Column(name='raw_Ha_med', format='E', array=np.array(rawHameds))
col33 = fits.Column(name='raw_r_med', format='E', array=np.array(rawrmeds))
col34 = fits.Column(name='moon_altitude', format='E', array=np.array(moon_altitudes))
col35 = fits.Column(name='moon_separation', format='E', array=np.array(moon_separations))
col36 = fits.Column(name='moon_phase', format='E', array=np.array(moon_phases))
col37 = fits.Column(name='r_exp_time', format='E', array=np.array(exptimesr))
cols = fits.ColDefs([col1, col2, col3, col4, col5, col37, col32, col33, col24, col22, col23, col6, col7, col8, col9, col10, col11, col12, col13, col14, col16, col17, col18, col19, col20, col21, col25, col26, col27, col28, col29, col30, col31, col34, col35, col36, col15])
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto(root + 'iphas-cleaning_'+str(Begin)+'.fits', clobber=True)

print ' '
print 'DONE.'
print ' '

