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

    ### Function make a cross section of a variable (or two!)
    ### Inputs:
    ###  var1, name of first variable to cross-section (contours)
    ###  var2, optional, name of second variable to cross-section (shading)
    ###  spath, optional, string, location to save the image
    ###  lonlat, location to take cross-section
    ###  direction, direction of cross-section. May be 'we' or 'sn'
    ###    for west-east and sotuh-north

    def cross_section(self, var1, lonlat, direction, var2=None, spath=None, xlimit=(-135, -60), ylimit=(1000, 50), nlevs=15):

        # Check if potetnial temp is desired
        theta_flag = False
        if (var1 == 'THETA'):
            var1 = 'TMP'
            theta_flag = True

        # Download the data
        plevs = [
            "40 mb", "50 mb", "70 mb", "100 mb", "150 mb", "200 mb", "250 mb", "300 mb", "350 mb",
            "400 mb", "450 mb", "500 mb", "550 mb", "600 mb", "650 mb", "700 mb", "750 mb", "800 mb",
            "850 mb", "900 mb", "925 mb", "950 mb", "975 mb", "1000 mb"
        ]
        if (var2 != None):
            self.set_vars({var1:plevs, var2: plevs, 'PRES': ['surface']})
        else:
            self.set_vars({var1:plevs, 'PRES': ['surface']})
        data = self.download_vars()

        # Locate the index of the grid for the cross-section
        # and extract data
        var1_cross = []
        # For west-east cross-sections
        if (direction.lower() == 'we'):
            ind = np.argmin((data['lats'][:,0] - lonlat)**2)
            x = data['lons'][ind,:]
            terrain = data['PRESsurface'][ind,:]/100.0 # hPa
            for p in plevs:
                var1_cross.append(data[var1+p][ind,:])

            # modify extent if necessary
            xlimit = list(xlimit)
            if (xlimit[0] < 0):
                xlimit[0] = xlimit[0]+360
            if (xlimit[1] < 0):
                xlimit[1] = xlimit[1]+360

            # Repeat for second variable if desired
            if (var2 != None):
                var2_cross = []
                for p in plevs:
                    var2_cross.append(data[var2+p][ind,:])

        # For south-north cross sections
        elif (direction.lower() == 'sn'):
            ind = np.argmin((data['lons'][:,0] - lonlat)**2)
            x = data['lats'][:,ind]
            terrain = data['PRESsurface'][ind,:]/100.0 # hPa
            for p in plevs:
                var1_cross.append(data[var1+p][:,ind])

            # Repeat for second variable if desired
            if (var2 != None):
                var2_cross = []
                for p in plevs:
                    var2_cross.append(data[var2+p][ind,:])

        else:
            raise ValueError("Direction must be either 'we' or 'sn'")

        # Convert to numpy arrays
        var1_cross = np.array(var1_cross)
        plevs = np.array([p.replace(' mb','') for p in plevs], dtype=np.float32)
        if (var2 != None):
            var2_cross = np.array(var2_cross)

        # Compute potential temperature if desired
        if theta_flag:
            var1_cross = var1_cross*(1000.0/np.stack([plevs]*x.size, axis=-1))**(287.0/1004.0)
            var1 = 'THETA (K)'

        # Nan the data outside of the bounds
        # X bounds
        mask = (np.stack([x]*plevs.size, axis=0)<xlimit[0]-0.30) | (np.stack([x]*plevs.size, axis=0)>xlimit[1]+0.30)
        var1_cross[mask] = np.nan
        if (var2 != None):
            var2_cross[mask] = np.nan

        # Y bounds
        mask = (np.stack([plevs]*x.size, axis=-1)>ylimit[0]+50) | (np.stack([plevs]*x.size, axis=-1)<ylimit[1]-50)
        var1_cross[mask] = np.nan
        if (var2 != None):
            var2_cross[mask] = np.nan

        # Make the figure
        fig, ax = pp.subplots(constrained_layout=True)

        cont1 = ax.contour(x, plevs, var1_cross, colors='red', levels=nlevs, zorder=9)
        ax.contour(x, plevs, var1_cross, colors='black', levels=nlevs, zorder=10)
        ax.clabel(cont1, inline=True, inline_spacing=0.5, zorder=11)

        # Shade terrain
        ax.fill_between(x, 1000, y2=terrain, color='sienna')

        # Second variable if desired
        if (var2 != None):
            cont2 = ax.contourf(x, plevs, var2_cross)
            cb = fig.colorbar(cont2, ax=ax, orientation='vertical')
            cb.set_label(var2)

        # Set Limits
        ax.set_xlim(xlimit)
        ax.set_ylim(ylimit)

        # Decorate
        ax.set_title(var1, fontsize=14, fontweight='bold')

        if (type(spath) == str):
                pp.savefig(spath)

        return fig, ax

    ### Function to create a QG PV map
    ### If save_data is a string, then the data will be saved as a numpy zip array to that path
    ### If load_file is a string, then that numpy zip file will be loaded.
    def qgpv_plot(self, spath=None, save_data=None, load_file=None, latlon_extent=[-135, -60, 20, 55], point=None, nbarbs=15, projection=ccrs.PlateCarree()):

        # Preserve the current variable list
        old_vars = self.vars

        # Grab variables necessary for QG PV
        self.set_vars({
            'ABSV': ['250 mb'],
            'HGT': ['250 mb'],
            'TMP': ['200 mb', '250 mb', '300 mb'],
        })

        # Retreive the data
        if (load_file == None):
            data = self.download_vars()
            
            # Check if saving data
            if (save_data != None):
                np.savez(save_data, **data)

        else:
            data = np.load(load_file)

        # Compute planetary vorticity
        pvort = 2.0*2.0*np.pi/(24.0*3600.0)*np.sin(data['lats']*np.pi/180.0)
        rvort = data['ABSV250 mb']-pvort

        # Compute the stability parameter
        theta_200mb = data['TMP200 mb']*((1000.0/200.0)**(287.05/1005.0))
        theta_300mb = data['TMP300 mb']*((1000.0/300.0)**(287.05/1005.0))
        Sp = -data['TMP250 mb'] * (np.log(theta_200mb)-np.log(theta_300mb))/(-100.0)

        # Compute thermal (i.e. stretching) vorticity
        tvort = (-(10**(-4))/Sp)*((np.log(data['TMP200 mb'])-np.log(data['TMP300 mb']))/(-100.0))

        # Compute QG PV
        qgpv = rvort+pvort+tvort

        # Mask out data outside region of interest
        mask = ~((data['lons']>latlon_extent[0]-0.1) & (data['lons']<latlon_extent[1]+0.1) & (data['lats']>latlon_extent[2]-0.1) & (data['lats']<latlon_extent[3]+0.1))
        #qgpv[mask] = np.nan
        #pvort[mask] = np.nan
        #rvort[mask] = np.nan
        #tvort[mask] = np.nan

        # Get boundaries of the data for consistent colormaps
        #bmin = min(np.nanmin(qgpv), np.nanmin(rvort), np.nanmin(pvort), np.nanmin(tvort))
        #bmax = max(np.nanmax(qgpv), np.nanmax(rvort), np.nanmax(pvort), np.nanmax(tvort))
        bmin = -0.001
        bmax = 0.001
        norm = mcolor.SymLogNorm(0.000000001, vmin=bmin, vmax=bmax)

        # Reset variable list
        self.vars = old_vars

        # Make a four panel with the total QG PV and each term
        fig, axes = pp.subplots(nrows=2, ncols=2, figsize=(18,12), subplot_kw={'projection': ccrs.PlateCarree()}, constrained_layout=True)
        fig.suptitle(f'Valid Time: {(self.date+timedelta(hours=self.forecast)).strftime("%Y-%m-%d %H:%M UTC")}', fontsize=16, fontweight='bold')

        # QG PV plot
        axes[0,0].set_title('QG PV - 250 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        cont = axes[0,0].pcolormesh(data['lons'], data['lats'], qgpv, norm=norm, cmap='coolwarm')

        # Thermal Vorticity
        axes[1,0].set_title('Thermal Vorticity', loc='left', ha='left', fontsize=16, fontweight='bold')
        axes[1,0].pcolormesh(data['lons'], data['lats'], tvort, norm=norm, cmap='coolwarm')

        # Relative Vorticity
        axes[0,1].set_title('Relative Vorticity', loc='left', ha='left', fontsize=16, fontweight='bold')
        axes[0,1].pcolormesh(data['lons'], data['lats'], rvort, norm=norm, cmap='coolwarm')

        # Planetary vorticity
        axes[1,1].set_title('Planetary Vorticity', loc='left', ha='left', fontsize=16, fontweight='bold')
        axes[1,1].pcolormesh(data['lons'], data['lats'], pvort, norm=norm, cmap='coolwarm')

        # Add a colorbar
        cb = fig.colorbar(cont, ax=list(axes.ravel()), orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Vorticity (s $^{-1}$)', fontsize=16, fontweight='bold')

        # Add decorations
        for ax in axes.flatten():

            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.BORDERS)
            ax.add_feature(cfeature.STATES)
            ax.set_extent(latlon_extent)
            ax.contour(data['lons'], data['lats'], data['HGT250 mb'], colors='black', levels=15)

        # Add a point of interest
        if (point != None):
            for ax in axes.flatten():
                ax.scatter(point[0], point[1], marker='x', s=30, c='darkgoldenrod', transform=ccrs.PlateCarree())

        # Save the plot if desired
        if (type(spath) == str):
            pp.savefig(spath)

        return

    ### Function to make a standardized 4-panel plot
    def four_panel(self, spath=None, latlon_extent=[-135, -60, 20, 55], point=None, nbarbs=15, projection=ccrs.PlateCarree()):

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
        data_proj = ccrs.PlateCarree()
        proj = projection

        # Make the plots
        fig, axes = pp.subplots(nrows=2, ncols=2, figsize=(18,12), subplot_kw={'projection':proj}, constrained_layout=True)

        fig.suptitle(f'Valid Time: {(self.date+timedelta(hours=self.forecast)).strftime("%Y-%m-%d %H:%M UTC")}', fontsize=16, fontweight='bold')

        # 250 mb
        axes[0,0].set_title('250 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        hcont = axes[0,0].contour(data['lons'], data['lats'], data['HGT250 mb'], colors='black', levels=20, transform=data_proj)
        axes[0,0].clabel(hcont, hcont.levels, inline=True, fontsize=10)

        norm = mcolor.BoundaryNorm(np.arange(30,125,5),self.wspd250.N, extend='both')
        fcont = axes[0,0].pcolormesh(data['lons'], data['lats'], wspd250, cmap=self.wspd250, norm=norm, shading='nearest', transform=data_proj)
        cb = fig.colorbar(fcont, ax=axes[0,0], orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Wind Speed (kts)', fontsize=14, fontweight='roman')

        axes[0,0].barbs(data['lons'][::nbarbs, ::nbarbs], data['lats'][::nbarbs, ::nbarbs],
                        data['UGRD250 mb'][::nbarbs,::nbarbs], data['VGRD250 mb'][::nbarbs,::nbarbs], transform=data_proj)

        # 500 mb
        axes[0,1].set_title('500 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        hcont = axes[0,1].contour(data['lons'], data['lats'], data['HGT500 mb'], colors='black', levels=20, transform=data_proj)
        axes[0,1].clabel(hcont, hcont.levels, inline=True, fontsize=10)

        norm = mcolor.BoundaryNorm(np.linspace(-0.0016, 0.0016, 33), self.vort500.N, extend='both')
        fcont = axes[0,1].pcolormesh(data['lons'], data['lats'], data['ABSV500 mb'], cmap=self.vort500, norm=norm, shading='nearest', transform=data_proj)
        cb = fig.colorbar(fcont, ax=axes[0,1], orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Absolute Vorticity (s$^{-1}$)', fontsize=14, fontweight='roman')

        axes[0,1].barbs(data['lons'][::nbarbs, ::nbarbs], data['lats'][::nbarbs, ::nbarbs],
                        data['UGRD500 mb'][::nbarbs,::nbarbs], data['VGRD500 mb'][::nbarbs,::nbarbs], transform=data_proj)

        # 700 mb
        axes[1,0].set_title('700 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        hcont = axes[1,0].contour(data['lons'], data['lats'], data['HGT700 mb'], colors='black', levels=20, transform=data_proj)
        axes[1,0].clabel(hcont, hcont.levels, inline=True, fontsize=10)

        norm = mcolor.BoundaryNorm(np.linspace(-20, 20, 21), self.temp700.N, extend='both')
        fcont = axes[1,0].pcolormesh(data['lons'], data['lats'], data['TMP700 mb']-273.15, cmap=self.temp700, norm=norm, shading='nearest', transform=data_proj)
        cb = fig.colorbar(fcont, ax=axes[1,0], orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Temperature (°C)', fontsize=14, fontweight='roman')

        mcont = axes[1,0].contour(data['lons'], data['lats'], data['RH700 mb'], colors='forestgreen', levels=[50,80], linestyles='--', transform=data_proj)
        axes[1,0].clabel(mcont, mcont.levels, inline=True, fontsize=10)

        axes[1,0].barbs(data['lons'][::nbarbs, ::nbarbs], data['lats'][::nbarbs, ::nbarbs],
                        data['UGRD700 mb'][::nbarbs,::nbarbs], data['VGRD700 mb'][::nbarbs,::nbarbs], transform=data_proj)

        # 850 mb
        axes[1,1].set_title('850 mb', loc='left', ha='left', fontsize=16, fontweight='bold')
        hcont = axes[1,1].contour(data['lons'], data['lats'], data['HGT850 mb'], colors='black', levels=20, transform=data_proj)
        axes[1,1].clabel(hcont, hcont.levels, inline=True, fontsize=10)

        norm = mcolor.BoundaryNorm(np.linspace(-20, 20, 21), self.temp700.N, extend='both')
        fcont = axes[1,1].pcolormesh(data['lons'], data['lats'], data['TMP850 mb']-273.15, cmap=self.temp700, norm=norm, shading='nearest', transform=data_proj)
        cb = fig.colorbar(fcont, ax=axes[1,1], orientation='horizontal', pad=0.02, aspect=50)
        cb.set_label('Temperature (°C)', fontsize=14, fontweight='roman')

        mcont = axes[1,1].contour(data['lons'], data['lats'], data['RH850 mb'], colors='forestgreen', levels=[50,80], linestyles='--', transform=data_proj)
        axes[1,1].clabel(mcont, mcont.levels, inline=True, fontsize=10)

        axes[1,1].barbs(data['lons'][::nbarbs, ::nbarbs], data['lats'][::nbarbs, ::nbarbs],
                        data['UGRD850 mb'][::nbarbs,::nbarbs], data['VGRD850 mb'][::nbarbs,::nbarbs], transform=data_proj)

        # Add decorations
        for ax in axes.flatten():

            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.BORDERS)
            ax.add_feature(cfeature.STATES)
            ax.set_extent(latlon_extent, crs=data_proj)

        # Add a point of interest
        if (point != None):
            for ax in axes.flatten():
                ax.scatter(point[0], point[1], marker='x', s=30, c='darkgoldenrod', transform=ccrs.PlateCarree())

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
