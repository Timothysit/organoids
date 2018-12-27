# get image dimensions of bar plots created before (to match them with the new figures)
from PIL import Image

im = Image.open('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/figures/paper_figures_first_draft/supplementary_figure/fireRate_bar_0413_slice4_withErrorBar_v2.eps')
width, height = im.size
dpi = 300
print(width / dpi)
print(height / dpi)