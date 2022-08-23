import glob
import os

gif_name = 'cross-rings'
file_list = glob.glob('gnomview*.png')
list.sort(file_list)

with open('image_list.txt', 'w') as file:
    for item in file_list:
        file.write("%s\n" % item)

os.system('convert @image_list.txt {}.gif'.format(gif_name))
