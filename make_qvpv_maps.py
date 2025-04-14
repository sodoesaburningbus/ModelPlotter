### This script creates QV PV maps
### Christopher Phillips

# Import modules
from plotter import model_plotter
from datetime import datetime
import cartopy.crs as ccrs

# The date
date = datetime(year=2024, month=5, day=3, hour=12)

# Create the object
model = model_plotter(date=date, source='GFS')

# Make the 4 panel plot
spath = f'qgpv_plot_{date.strftime("%Y%m%d_%H%MUTC")}.png'
model.qgpv_plot(spath=spath)#, load_file='qgpv_test.npz')