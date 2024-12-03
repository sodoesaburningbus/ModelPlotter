### This script creates the necessary plots for the MET 480 Lab
### Christopher Phillips

# Import modules
from plotter import model_plotter
from datetime import datetime
import cartopy.crs as ccrs

# The date
date = datetime(year=2021, month=1, day=31, hour=12)

# Create the object
model = model_plotter(date=date, source='GFS')

# Make the 4 panel plot
spath = f'cross_section_{date.strftime("%Y%m%d_%H%MUTC")}.png'
model.cross_section('THETA', 40, 'we', spath=spath, nlevs=20, ylimit=(1000,200))