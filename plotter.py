### This object contains methods for creating
### plots using common weather model data sources
### such as HRRR and GFS.
###
### Christopher Phillips
### Valparaiso University
### October 2024

### Import required modules
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta
import matplotlib.pyplot as pp
import matplotlib.colors as mcolor
import numpy as np
import os
import pygrib
import urllib.request as urlreq

### The main object
### Inputs:
###  source, optional, string, the model to use.
###  date, optional, datetime object, initialization date and time of the desired model run (UTC)
###  vars, optional, dictionary with variable names as keys and desired levels as a list of strings
###    e.g. {'TMP': ['1000 mb','850 mb']}
###  forecast, optional, integer, the forecast hour to use (date is the initialization time, default=0)
class model_plotter:

    ### Initialize the function
    def __init__(self, source='GFS', date=None, vars=[], forecast=0):

        # Declare the base URLs
        self.base_urls = {
            'gfs': 'https://noaa-gfs-bdp-pds.s3.amazonaws.com/',
            'hrrr': 'https://noaa-hrrr-bdp-pds.s3.amazonaws.com/'
        }

        # Set attributes with the input data
        self.source = source.lower()
        self.date = date
        self.forecast = forecast
        self.vars = vars

        # Update file URLs
        self.update_urls()

        # Create the colormaps for standadized plots
        self.create_colormaps()

        return
    
    ### Function to download the requested variables
    def download_vars(self):

        # Create a dictionary to hold the data
        data = {}

        # loop over each variable and level
        for var in self.vars.keys():

            for lev in self.vars[var]:

                # Check the index
                flag = False
                for line in self.index:

                    # Check if var already found
                    if flag:
                        byte_end = line.split(":")[1]
                        break

                    # Locate variable
                    if ((var in line) and (lev in line)):
                        byte_start = line.split(":")[1]
                        flag = True

                # Download the data
                os.system("curl -o model_out -r {}-{} {}".format(byte_start, byte_end, self.url+self.file))
                grib = pygrib.open("model_out")
                data[var+lev] = grib[1].values
                lons = grib[1].longitudes.reshape(data[var+lev].shape)
                lats = grib[1].latitudes.reshape(data[var+lev].shape)
                grib.close()

        # Fix any negative longitudes
        lons = np.where(lons<0, lons+360, lons)

        # Save the latitudes and longitudes
        data['lons'] = lons
        data['lats'] = lats

        return data

    ### Function to make a single plot of each variable


    ### Function to make a standardized 4-panel plot
    def four_panel(self, spath=None, latlon_extent=[-135, -60, 20, 55], point=None):

        # Preserve the current variable list
        old_vars = self.vars

        # Set the variables to the standard list
        self.set_vars({
            'TMP': ['850 mb', '700 mb'],
            'HGT': ['850 mb','700 mb','500 mb','250 mb'],
            'RH': ['850 mb', '700 mb'],
            'UGRD': ['850 mb','700 mb','500 mb','250 mb'],
            'VGRD': ['850 mb','700 mb','500 mb','250 mb'],
            'ABSV': ['500 mb']
        })

        # Download the data
        data = self.download_vars()

        # Compute wind speed at 250 mb
        wspd250 = np.sqrt(data['UGRD250 mb']**2+data['VGRD250 mb']**2)*1.94

        # Create the map projection
        proj = ccrs.PlateCarree()

        # Make the plots
        fig, axes = pp.subplots(nrows=2, ncols=2, figsize=(18,12), subplot_kw={'projection':proj}, constrained_layout=True)

        fig.suptitle(f'Valid Time: {(self.date+timedelta(hours=self.forecast)).strftime("%Y-%m-%d %H:%M UTC")}', fontsize=16, fontweight='bold')

        # 250 mb
        axes[0,0].set_title('250 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        hcont = axes[0,0].contour(data['lons'], data['lats'], data['HGT250 mb'], colors='black', levels=20)
        axes[0,0].clabel(hcont, hcont.levels, inline=True, fontsize=10)

        norm = mcolor.BoundaryNorm(np.arange(30,125,5),self.wspd250.N, extend='both')
        fcont = axes[0,0].pcolormesh(data['lons'], data['lats'], wspd250, cmap=self.wspd250, norm=norm, shading='nearest')
        cb = fig.colorbar(fcont, ax=axes[0,0], orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Wind Speed (kts)', fontsize=14, fontweight='roman')

        # 500 mb
        axes[0,1].set_title('500 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        hcont = axes[0,1].contour(data['lons'], data['lats'], data['HGT500 mb'], colors='black', levels=20)
        axes[0,1].clabel(hcont, hcont.levels, inline=True, fontsize=10)

        norm = mcolor.BoundaryNorm(np.linspace(-0.0016, 0.0016, 33), self.vort500.N, extend='both')
        fcont = axes[0,1].pcolormesh(data['lons'], data['lats'], data['ABSV500 mb'], cmap=self.vort500, norm=norm, shading='nearest')
        cb = fig.colorbar(fcont, ax=axes[0,1], orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Absolute Vorticity (s$^{-1}$)', fontsize=14, fontweight='roman')

        # 700 mb
        axes[1,0].set_title('700 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        hcont = axes[1,0].contour(data['lons'], data['lats'], data['HGT700 mb'], colors='black', levels=20)
        axes[1,0].clabel(hcont, hcont.levels, inline=True, fontsize=10)

        norm = mcolor.BoundaryNorm(np.linspace(-20, 20, 21), self.temp700.N, extend='both')
        fcont = axes[1,0].pcolormesh(data['lons'], data['lats'], data['TMP700 mb']-273.15, cmap=self.temp700, norm=norm, shading='nearest')
        cb = fig.colorbar(fcont, ax=axes[1,0], orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Temperature (°C)', fontsize=14, fontweight='roman')

        mcont = axes[1,0].contour(data['lons'], data['lats'], data['RH700 mb'], colors='forestgreen', levels=[50,80], linestyles='--')
        axes[1,0].clabel(mcont, mcont.levels, inline=True, fontsize=10)

        # 850 mb
        axes[1,1].set_title('850 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        hcont = axes[1,1].contour(data['lons'], data['lats'], data['HGT850 mb'], colors='black', levels=20)
        axes[1,1].clabel(hcont, hcont.levels, inline=True, fontsize=10)

        norm = mcolor.BoundaryNorm(np.linspace(-20, 20, 21), self.temp700.N, extend='both')
        fcont = axes[1,1].pcolormesh(data['lons'], data['lats'], data['TMP850 mb']-273.15, cmap=self.temp700, norm=norm, shading='nearest')
        cb = fig.colorbar(fcont, ax=axes[1,1], orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Temperature (°C)', fontsize=14, fontweight='roman')

        mcont = axes[1,1].contour(data['lons'], data['lats'], data['RH850 mb'], colors='forestgreen', levels=[50,80], linestyles='--')
        axes[1,1].clabel(mcont, mcont.levels, inline=True, fontsize=10)

        # Add decorations
        for ax in axes.flatten():

            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.BORDERS)
            ax.add_feature(cfeature.STATES)
            ax.set_extent(latlon_extent)

        # Add a point of interest
        if (point != None):
            for ax in axes.flatten():
                ax.scatter(point[0], point[1], marker='x', s=30, c='black', transform=ccrs.PlateCarree())

        # Save the plot if desired
        if (type(spath) == str):
            pp.savefig(spath)

        # Reset variable list
        self.vars = old_vars

        return fig, axes


    ### Function to make a standardized surface analysis plot

    ### Function to create default colormaps
    def create_colormaps(self):

        # 850 colormaps
        colors = {
            0.0: "#8B00FF",   # Violet (at -35)
            0.25: "#0000FF",  # Blue (at -20)
            0.5833: "#008080",  # Teal (at 0)
            0.75: "#00FF00",  # Green (at 10)
            0.9167: "#FFFF00",  # Yellow (at 20)
            1.0: "#FF0000"    # Red (at 25)
        }

        self.temp850 = mcolor.LinearSegmentedColormap.from_list('temp850', list(colors.items()))

        # 700 colormaps
        colors = {
            0.0: "#8B00FF",   # Violet (at -20)
            0.25: "#0000FF",  # Blue (at -10)
            0.5: "#008080",   # Teal (at 0)
            0.8: "#FFFF00",   # Yellow (at 12)
            1.0: "#FF0000"    # Red (at 20)
        }
        self.temp700 = mcolor.LinearSegmentedColormap.from_list('temp700', list(colors.items()))

        # 500 colormaps
        self.vort500 = pp.get_cmap('RdGy_r')

        # 250 colormap
        colors = {
            0.0: "#FFFFFF",   # White (at 30)
            0.25: "#AAFFFF",  # Light Teal (intermediate)
            0.5: "#008080",   # Teal
            0.7: "#9370DB",   # Lavender (transition from teal to pink)
            0.85: "#FF69B4",  # Pink
            1.0: "#C71585"    # Deeper Pink (at 120)
        }
        self.wspd250 = mcolor.LinearSegmentedColormap.from_list('wspd250', list(colors.items()))

        return

    ### Function to change the date
    ### Input
    ###   date, datetime object, initialization date and time of the desired model run (UTC)
    ###   forecast, optional, integer, the forecast hour to use (default=0)
    def set_date(self, date, forecast=0):

        # Change time info
        self.date = date
        self.forecast = forecast

        # Update the URLs
        self.update_urls()

        return

    ### Function to reset the variable list
    def set_vars(self, vars):

        self.vars = vars
        return

    ### Function create URLs
    ### This function also pulls in the index file for the data source
    def update_urls(self):

        # update URLs
        if (self.source == 'gfs'):
            self.url = f'{self.base_urls[self.source]}gfs.{self.date.strftime("%Y%m%d")}/{self.date.hour:02d}/'
            self.index_file = f'gfs.t{self.date.hour:02d}z.pgrb2.0p25.f{self.forecast:03d}.idx'
            self.file = f'gfs.t{self.date.hour:02d}z.pgrb2.0p25.f{self.forecast:03d}'

        # Read the model index file from the index url
        print(f'Downloading data from: {self.url+self.index_file}')
        self.index = []
        try:
            for line in urlreq.urlopen(self.url+self.index_file):
                self.index.append(line.decode('utf-8'))
        except:
            self.index = []
            if (self.source == 'gfs'):
                self.url = self.url+'atmos/'
                for line in urlreq.urlopen(self.url+self.index_file):
                    self.index.append(line.decode('utf-8'))

        return
