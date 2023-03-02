import sys
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as apc
import astropy.units as apu
import pandas as pd
import sqlalchemy

sys.path.append('/Users/amsys/paper_tss')
import tfpy

db = tfpy.open_database('qso_100', "root", "password")

def query_line(line):
    query = session.query(
        tfpy.Photon.Wavelength, tfpy.Photon.Resonance, tfpy.Photon.Delay, tfpy.Photon.Weight
    ).filter(tfpy.Photon.Resonance == line)
    data = pd.DataFrame(
        query.all(), columns=['Wavelength', 'Resonance', 'Delay', 'Weight']
    )
    return np.average(data['Wavelength'], weights=data['Weight'])

Session = sqlalchemy.orm.sessionmaker(bind=db)
session = Session()
query = session.query(tfpy.Photon.Delay, tfpy.Photon.Weight).filter(tfpy.Photon.Resonance.in_([131, 132]))
data = pd.DataFrame(
    query.all(), columns=['Delay', 'Weight']
)
delay = (np.average(data['Delay'], weights=data['Weight']) * apu.second)

tf = tfpy.TransferFunction(db, 'qso_100_mg', continuum=1,
            wave_bins=200, delay_bins=200)
tf.lines([131, 132]).run()
delay_2 = tf.delay(days=False)

print(delay, delay_2)
# Mg II: 2799.117
# Lines between 2798-2800: 130 (few), 131, 132, 139 (few)
# All 131-132 are Mg II

# C III: 1908.734
# Lines between 1908-1909: 323, 324, 325, 326, 200130
#
