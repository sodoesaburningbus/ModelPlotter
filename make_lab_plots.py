### This script creates the necessary plots for the MET 480 Lab
### Christopher Phillips

# Import modules
from plotter import model_plotter
from datetime import datetime

# The date
date = datetime(year=2022, month=3, day=30, hour=12)

# Create the object
model = model_plotter(date=date, source='GFS')

# Make the 4 panel plot
model.four_panel(spath=f'four_panel_{date.strftime("%Y%m%d_%H%MUTC")}.png', point=(-86.81, 33.52))