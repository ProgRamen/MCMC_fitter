#This piece of code takes a bunch of images in a directory
#then bins each image by a specified factor, before
#placing the resultant binned image in another chosen directory.
#This code is used to bin the cleaned images before mosaicking
#in the cleaning code.
import os


#Root folder containg input and output directories 
#as well as the factor by which to bin by are
#required inputs when calling the code.
def binning(root_folder, bin_factor):
	'''bin multiple files using montage'''


#Specify the input directory containing the images to bin
	input_dir = root_folder + 'Ha_cleandir/'
#Specify the output directory to place the binned imaged
	output_dir = root_folder + 'rawdir/'

#Uses montage command mShrink to bin.
	for filename in os.listdir(input_dir):
		os.system("mShrink " + input_dir + filename + " " + output_dir + filename + " " + str(bin_factor) + " > /dev/null 2>&1 &")