import matplotlib.pyplot as plt
import pandas as pd

from astropy.coordinates import SkyCoord, EarthLocation, Angle
from astropy.time import Time
import astropy.units as u
from astroplan import Observer, FixedTarget, is_observable, observability_table
from astroplan import AtNightConstraint, AltitudeConstraint
from astroplan.plots import plot_airmass


def define_observer():
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


def targetsFromCSV(path):
    """
    Obtain a target list from a csv file.

    Reads a csv file with the TESS alerts as downloaded from the TEV homepage,
    extracts coordinates and names, and transforms the entries to astroplan
    target objects.

    Parameters
    ----------
    path : str
        path to the csv file

    Returns
    -------
    alerts : pandas DataFrame
        table containing alerts
    targets : list
        list containing astroplan targets
    """
    alerts = pd.read_csv(path)
    targets = []
    for toi in alerts.iterrows():
        target = toi[1]
        coords = SkyCoord(target.RA, target.Dec, unit=(u.deg, u.deg))
        targets.append(FixedTarget(coord=coords, name=target['toi_id']))
    return alerts, targets


def define_constraints(minAltitude, start_time, end_time):
    """
    Define various constraints to define 'observability'.

    Format for times: YYYY-MM-DD

    Parameters
    ----------
    minAltitude : int
        minimum local elevation in degree
    start_time : str
        Earliest time of observation
    end_time : str
        Latest time of observation

    Returns
    -------
    constraints : list
        list of constraints
    start_time : Time object
        Earliest time of observation
    end_time : Time object
        Latest time of observation
    """
    constraints = [AtNightConstraint.twilight_astronomical(),
                   AltitudeConstraint(min=minAltitude*u.deg)]
    return constraints, Time(start_time), Time(end_time)


def check_observability(constraints, observer, targets, start_time, end_time):
    """
    Check observability from observatory given the observational constraints.

    Parameters
    ----------
    constraints : list
        list of constraints
    observer : astroplan Observer object
        e.g. Calar Alto Observatory, Spain
    targets : list
        list containing astroplan targets
    start_time : Time object
        Earliest time of observation
    end_time : Time object
        Latest time of observation

    Returns
    -------
    alerts : pandas DataFrame
        filtered alerts table containing observable targets
    """
    observabilityMask = is_observable(constraints, observer=CAHA,
                                  targets=targets, time_range=[start_time, end_time])
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

def defaultRun():
    CAHA = define_observer()
    alerts, targets = targetsFromCSV('data/toi-2019-01-25.csv')
    constraints, start_time, end_time = define_constraints(30,
        '2019-02-08 12:00', '2019-12-31 12:00')
    alerts = check_observability(constraints, CAHA, targets, start_time,
                                 end_time)


if __name__ == "__main__":
    defaultRun()
