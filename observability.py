import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from astropy.coordinates import SkyCoord, EarthLocation, Angle
from astropy.time import Time
import astropy.units as u
from astroplan import Observer, FixedTarget, is_observable, observability_table
from astroplan import AtNightConstraint, AltitudeConstraint
from astroplan.plots import plot_airmass
from astroplan.utils import time_grid_from_range


def defineCAHA():
    """
    Define Calar Alto observational site.

    Parameters
    ----------
    None

    Returns
    -------
    CAHA : astroplan Observer object
        Calar Alto Observatory, Spain
    """
    return Observer(longitude=Angle('-2d32.8m'), latitude=Angle('37d13.4m'),
                    elevation=2168*u.m, name="CAHA", timezone="Europe/Madrid")


def targetsFromCSV(path, minUpdated=None, **read_csv_kwargs):
    """
    Obtain a target list from a csv file.

    Reads a csv file with the TESS alerts as downloaded from the TEV homepage,
    extracts coordinates and names, and transforms the entries to astroplan
    target objects.

    Parameters
    ----------
    path : str
        path to the csv file
    minUpdated : str, optional
        earliest time in column "Updated" to be considered
        Format: YYYY-MM-DD
    read_csv_kwargs : keyword arguments, optional
        keyword arguments to pass to pandas' read_csv function

    Returns
    -------
    alerts : pandas DataFrame
        table containing alerts
    targets : list
        list containing astroplan targets
    """
    alerts = pd.read_csv(path, **read_csv_kwargs)
    alerts.rename(columns={'Full TOI ID': 'TOI'}, inplace=True)
    targets = []

    if minUpdated is not None:
        if not 'Updated' in alerts.columns:
            # newer alert files contain "Alerted" instead of "Updated"
            alerts.rename(columns={'Alerted':'Updated'}, inplace=True)

        alerts.Updated = pd.to_datetime(alerts.Updated)
        alerts.Updated = alerts.Updated.dt.tz_localize(None) # strip timezone
        minUpdated = np.datetime64(minUpdated)
        alerts = alerts[alerts.Updated > minUpdated]

    for toi in alerts.iterrows():
        target = toi[1]
        coords = SkyCoord(target['TIC Right Ascension'], target['TIC Declination'],
                          unit=(u.deg, u.deg))
        targets.append(FixedTarget(coord=coords, name=target['TOI']))
    return alerts, targets


def define_constraints(minAltitude, earliestObs, latestObs):
    """
    Define various constraints to define 'observability'.

    Format for times: YYYY-MM-DD

    Parameters
    ----------
    minAltitude : int
        minimum local elevation in degree
    earliestObs : str
        Earliest time of observation
    latestObs : str
        Latest time of observation

    Returns
    -------
    constraints : list
        list of constraints
    earliestObs : Time object
        Earliest time of observation
    latestObs : Time object
        Latest time of observation
    """
    constraints = [AtNightConstraint.twilight_astronomical(),
                   AltitudeConstraint(min=minAltitude*u.deg)]
    return constraints, Time(earliestObs), Time(latestObs)


def check_observability(alerts, constraints, observer, targets, earliestObs,
                        latestObs):
    """
    Check observability from observatory given the observational constraints.

    Parameters
    ----------
    alerts : pandas DataFrame
        table containing alerts
    constraints : list
        list of constraints
    observer : astroplan Observer object
        e.g. Calar Alto Observatory, Spain
    targets : list
        list containing astroplan targets
    earliestObs : Time object
        Earliest time of observation
    latestObs : Time object
        Latest time of observation

    Returns
    -------
    alerts : pandas DataFrame
        filtered alerts table containing observable targets
    """
    observabilityMask = is_observable(constraints, observer,
                                  targets, time_range=[earliestObs, latestObs])
    return alerts[observabilityMask]


def MdwarfFilter(candidates, maxTeff=3800):
    """
    Filter candidates list for M dwarfs.

    Parameters
    ----------
    candidates : pandas DataFrame
        Table containing candidates. Must contain a column "Teff".
    maxTeff : int
        maximum effective temperature in Kelvin

    Returns
    -------
    candidates : pandas DataFrame
        filtered candidates Table containing only M dwarfs
    """
    Mcandidates = candidates[candidates.Teff < maxTeff]


def plot_observability(candidates, constraints, observer, earliestObs,
                       latestObs, timeRes=24*u.hour, timeSubRes=.2*u.hour,
                       fig=None, ax=None, **kwargs):
    """ Visualize long-term observability of a set of targets.

    Parameters
    ----------
    candidates : pandas DataFrame
        Table containing candidates. Must contain a column "Teff".
    constraints : list
        list of constraints
    observer : astroplan Observer object
        e.g. Calar Alto Observatory, Spain
    earliestObs : Time object
        Earliest time of observation
    latestObs : Time object
        Latest time of observation
    timeRes : quantity
        time-resolution of the resulting plot, should be greater than 24h
    timeSubRes : quantity
        fine resolution used for computation of observable time fractions at
        each resolved time in the plot. Has to be (more than an order of
        magnitude) smaller than 'timeRes' for reasonable results.
    fig : matplotlib figure object, optional
        figure to plot on
    ax : matplotlib axis object, optional
        axis to plot on
    **kwargs : keyword arguments to pass to matplotlib's 'imshow'

    Returns
    -------
    fig : matplotlib figure
        figure containing the plot
    ax : matplotlib axis
        axis containing the plot
    """

    # flip candidate list
    candidates = candidates[::-1]

    # prepare the grid
    rough_grid = time_grid_from_range([earliestObs, latestObs],
                                     time_resolution=timeRes)
    observability_grid = np.zeros((len(candidates), len(rough_grid) - 1))
    for i, t in enumerate(rough_grid[:-1]):
        obsTable = observability_table(constraints, observer, candidates,
                                       time_grid_from_range([rough_grid[i],
                                       rough_grid[i+1]],
                                     time_resolution=timeSubRes))
        hoursObservable = obsTable['fraction of time observable']*24*u.hour
        observability_grid[:,i] = hoursObservable

    # now for the actual plot
    if ax == None:
        fig, ax = plt.subplots()
    dates = rough_grid.plot_date
    extent = [dates[0], dates[-1], -.5, len(candidates) -.5]
    im = ax.imshow(observability_grid, cmap='inferno', aspect='auto',
                   extent=extent, origin='lower', **kwargs)

    # get date ticks right
    ax.xaxis_date()
    date_format = mdates.DateFormatter('%b \'%y')
    ax.xaxis.set_major_formatter(date_format)
    ax.tick_params(rotation=45, axis='x', which='both', length=3, direction='out',
                   pad=2)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_horizontalalignment('left')

    # some eyecandy
    ax.set_yticks(range(len(candidates)))
    ax.set_yticklabels([c.name for c in candidates])
    ax.set_xlabel('Date', labelpad=9)
    ax.set_ylabel('TOI')
    ax.grid(False)
    fig.colorbar(im).set_label('Hours Observable per Night', labelpad=15)

    return fig, ax


def defaultRun():
    CAHA = CAHA()
    alerts, targets = targetsFromCSV('data/toi-2019-01-25.csv')
    constraints, earliestObs, latestObs = define_constraints(30,
        '2019-02-08 12:00', '2019-12-31 12:00')
    observables = check_observability(alerts, constraints, CAHA, targets, earliestObs,
                                 latestObs)


if __name__ == "__main__":
    defaultRun()
