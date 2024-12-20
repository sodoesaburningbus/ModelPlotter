### This script creates the necessary plots for the MET 480 Lab
### Christopher Phillips

# Import modules
from plotter import model_plotter
from datetime import datetime
import cartopy.crs as ccrs

# The date
date = datetime(year=2024, month=10, day=20, hour=12)

# Create the object
model = model_plotter(date=date, source='GFS')

# Make the 4 panel plot
spath = f'four_panel_{date.strftime("%Y%m%d_%H%MUTC")}.png'
model.four_panel(spath=spath, point=(-88.11, 44.5), latlon_extent=[-180, -80, 35, 85])