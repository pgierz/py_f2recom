"""
Visualization tools for FESOM2-REcoM2 model output.
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def create_basemap(projection=ccrs.PlateCarree(), figsize=(10, 6)):
    """
    Create a basic map for plotting model output.
    
    Parameters
    ----------
    projection : cartopy.crs.Projection, optional
        Map projection to use (default: PlateCarree)
    figsize : tuple, optional
        Figure size in inches (width, height)
        
    Returns
    -------
    tuple
        (figure, axis) objects
    """
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=projection)
    
    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    
    # Add gridlines
    if projection != ccrs.PlateCarree():
        ax.gridlines()
    
    return fig, ax

def plot_global_field(data, lon, lat, title=None, cmap='viridis', **kwargs):
    """
    Plot a global field on a map.
    
    Parameters
    ----------
    data : array-like
        2D array of data to plot
    lon : array-like
        Longitude coordinates
    lat : array-like
        Latitude coordinates
    title : str, optional
        Plot title
    cmap : str, optional
        Colormap to use
    **kwargs : dict, optional
        Additional arguments passed to pcolormesh
    """
    fig, ax = create_basemap()
    
    # Plot the data
    mesh = ax.pcolormesh(lon, lat, data, transform=ccrs.PlateCarree(),
                        cmap=cmap, **kwargs)
    
    # Add colorbar
    plt.colorbar(mesh, ax=ax, orientation='vertical',
                label=kwargs.get('label', ''))
    
    if title:
        ax.set_title(title)
    
    plt.tight_layout()
    return fig, ax
