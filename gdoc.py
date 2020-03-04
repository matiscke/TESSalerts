""" Script for interaction between a google spreadsheet and pandas dataframes.

adapted from https://gist.github.com/dcbark01/8f1240d2353b44b7f9a39871154bbe5f#file-google_sheet_tut-py
by Daniel Barker.
"""

from apiclient.discovery import build
from httplib2 import Http
from oauth2client import file, client, tools
import pandas as pd


def get_google_sheet(spreadsheet_id, range_name):
    """ Retrieve sheet data using OAuth credentials and Google Python API. """
    scopes = 'https://www.googleapis.com/auth/spreadsheets.readonly'
    # Setup the Sheets API
    store = file.Storage('credentials.json')
    creds = store.get()
    if not creds or creds.invalid:
        flow = client.flow_from_clientsecrets('client_secret.json', scopes)
        creds = tools.run_flow(flow, store)
    service = build('sheets', 'v4', http=creds.authorize(Http()))

    # Call the Sheets API
    gsheet = service.spreadsheets().values().get(spreadsheetId=spreadsheet_id, range=range_name).execute()
    return gsheet


def gsheet2df(gsheet):
    """ Converts Google sheet data to a Pandas DataFrame.
    Note: This script assumes that your data contains a header file on the THIRD row!
    Also note that the Google API returns 'none' from empty cells - in order for the code
    below to work, you'll need to make sure your sheet doesn't contain empty cells,
    or update the code to account for such instances.
    """
    header = gsheet.get('values', [])[2]   # Assumes THIRD line is header!
    values = gsheet.get('values', [])[3:]  # Everything else is data.
    if not values:
        print('No data found.')
    else:
        all_data = []
        for col_id, col_name in enumerate(header):
            column_data = []
            for row in values:
                try:
                    column_data.append(row[col_id])
                except IndexError:
                    break
            ds = pd.Series(data=column_data, name=col_name)
            all_data.append(ds)
        df = pd.concat(all_data, axis=1)
        return df


def read_sheet(sheet_ID, range_name='Candidates'):
    """ read spreadsheet into a pandas dataframe"""
    gsheet = get_google_sheet(sheet_ID, range_name)
    df = gsheet2df(gsheet)
    print('Dataframe size = ', df.shape)
    return df


