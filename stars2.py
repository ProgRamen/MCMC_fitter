#This code uses tha Tycho catalogue of bright stars
#to find and list all cases where bright stars fall on 
#or near IPHAs images between a pre-defined magnitude range.

#Required files in same directory as code: 
#'iphas-images.fits'
#'tycho____.fit' where _____ is whatever region you're looking at.

from astropy.io import fits
import numpy as np


#Specify the minimum and maximum magnitudes to search for.
magMin = raw_input('Min Mag: ')
magMax = raw_input('Max Mag: ')


#Calling the list of stars (usually tycho.fit). To search only a chosen area
#of sky, this list must be masked to include only those in
#the region desired before running this code.
catalogue = fits.open('tychoHEART.fit')[1].data		

#The next 4 lines call the desired collums from the table
#into arrays
mags = catalogue.field(3) #hipparcos3.fit - 3,0,1
raPos = catalogue.field(0) #tycho3.fit - 3,0,1
decPos = catalogue.field(1)
name = catalogue.field(2)
print len(mags)

#Calling the list of IPHAS images.
data = fits.open('iphas-images.fits')[1].data

#The next 8 lines call the required collumns from the IPHAS
#table into arrays
raMin = data.field(20)
raMax = data.field(21)
decMin = data.field(22)
decMax = data.field(23)
runs = data.field(0)
ccds = data.field(1)
fields = data.field(7)
urls = data.field(2)


#Next, the IPHAS images are masked to only incude the best
#data available for each pointing.
bests = data.field(9)
bests = bests == True
raMin = raMin[bests]
decMin = decMin[bests]
raMax = raMax[bests]
decMax = decMax[bests]
runs = runs[bests]
ccds = ccds[bests]
fields = fields[bests]
urls = urls[bests]




#Now the table of stars is masked to only include stars 
#within the chosen magnitude range.
raRange=[]
decRange=[]
magRange=[]
nameRange=[]
count=0
i=0 
while i<len(mags):
	mag = mags[i]
	raP = raPos[i]
	decP = decPos[i]
	nameP = name[i]
	if int(magMin)<=mag<int(magMax):
		magRange.append(mag)
		raRange.append(raP)
		decRange.append(decP)
		nameRange.append(nameP)
	i+=1
magRange = np.array(magRange)
raRange = np.array(raRange)
decRange = np.array(decRange)
nameRange = np.array(nameRange)



inImage = []
one = []
two = []
three = []
four = []
five = []
six = []
seven = []
eight = []
nine = []
ten = []
eleven = []




#The main body of the code.
save = 1000
starCount = 0
counting = 0
j=0
while j<len(magRange):
#The iteration is printed to the command line.
	print str(j+1)+'/'+str(len(magRange))

#This iterations star is called in
	ra = raRange[j]
	dec = decRange[j]
	magc = magRange[j]
	namec = nameRange[j]


#Some stars in the table are wrongly marked as having 0 mag.
#These are masked.
	if magc == 0.:
		zeromag = '*'
		magc = 0.01
	else:
		zeromag = ' '




	k=0
	while k<len(raMax):

#Each image is then checked for proximity to the star.
		ramin=raMin[k]
		ramax=raMax[k]
		decmin=decMin[k]
		decmax=decMax[k]
		run = runs[k]
		ccd = ccds[k]
		field = fields[k]
		url = urls[k]

#The radius of the star according to its magnitude is calculated.
#This is used to see if the star will have a affect on an image
#even if it only falls nearby.
		radius = int((8192/(float(magc)**2))+((1/float(magc))*1000)+100)

		radius_arcsec = float(radius)*0.33
		r = float(radius)*9.167e-5
		
#The conditions for a star or its surrounding scattered light 
#falling within the IPHAS image are set.
		cond1 = ra<(ramax+r)
		cond2 = ra>(ramin-r)
		cond3 = dec<(decmax+r)
		cond4 = dec>(decmin-r)
		inImageWhen = cond1&cond2&cond3&cond4


		if inImageWhen == True:
#If the star does fall on or near the image, all information
#is saved.
			count+=1
			one.append(str(run))
			two.append(str(ccd))
			inImage.append(str(run)+'-'+str(ccd))
			three.append(str(magc))
			four.append(str(ra))
			five.append(str(dec))
			six.append(str(namec))
			seven.append(str(field))
			eight.append(str(url))
			nine.append(radius)
			ten.append(radius_arcsec)
			eleven.append(str(zeromag))
			#print str(run)+'-'+str(ccd)
		else:
			count=count
		k+=1
#the counts keep track of the number of images affected, and
#how many stars affect them.
	if counting<count:
		starCount+=1
		counting = count

#The output table of stars in images is saves with '_PART' 
#at the end of its filename every 1000 iterations to keep track
#or progress.
	if int(j) == int(save):
		print 'saving progress...'
		col1 = fits.Column(name='run', format='20A', array=np.array(one))
		col2 = fits.Column(name='ccd', format='20A', array=np.array(two))
		col3 = fits.Column(name='star mag', format='20A', array=np.array(three))
		col4 = fits.Column(name='star ra', format='E', array=np.array(four))
		col5 = fits.Column(name='star dec', format='E', array=np.array(five))
		col6 = fits.Column(name='star name (TYC)', format='20A', array=np.array(six))
		col7 = fits.Column(name='field id', format='20A', array=np.array(seven))
		col8 = fits.Column(name='image url', format='55A', array=np.array(eight))
		col9 = fits.Column(name='r (pix)', format='E', array=np.array(nine))
		col10 = fits.Column(name='r (arcsec)', format='E', array=np.array(ten))
		col11 = fits.Column(name='no Vmag data', format='5A', array=np.array(eleven))
		cols = fits.ColDefs([col1, col2, col7, col6, col3, col4, col5, col8, col9, col10, col11])
		tbhdu = fits.BinTableHDU.from_columns(cols)
		tbhdu.writeto('HEARTStars'+str(magMin)+'-'+str(magMax)+'_PART.fits', clobber=True)
		save = save + 1000



	j+=1

#When all stars and images are tested, all cases where a star
#falls on or near an image are stored to a fits table.
col1 = fits.Column(name='run', format='20A', array=np.array(one))
col2 = fits.Column(name='ccd', format='20A', array=np.array(two))
col3 = fits.Column(name='star mag', format='20A', array=np.array(three))
col4 = fits.Column(name='star ra', format='E', array=np.array(four))
col5 = fits.Column(name='star dec', format='E', array=np.array(five))
col6 = fits.Column(name='star name (TYC)', format='20A', array=np.array(six))
col7 = fits.Column(name='field id', format='20A', array=np.array(seven))
col8 = fits.Column(name='image url', format='55A', array=np.array(eight))
col9 = fits.Column(name='r (pix)', format='E', array=np.array(nine))
col10 = fits.Column(name='r (arcsec)', format='E', array=np.array(ten))
col11 = fits.Column(name='no Vmag data', format='5A', array=np.array(eleven))
cols = fits.ColDefs([col1, col2, col7, col6, col3, col4, col5, col8, col9, col10, col11])
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto('HEARTStars'+str(magMin)+'-'+str(magMax)+'.fits', clobber=True)


#The final tally of the number of stars the fall on or near,
#images and the number of images affected is output to the
#command line.
print 'For magnitude range '+str(magMin)+'-'+str(magMax)+', there are '+str(starCount)+' stars in '+ str(count)+' images.'





