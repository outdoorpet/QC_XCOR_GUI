import pyasdf
from obspy.clients.fdsn.client import Client
from obspy.clients.fdsn.header import FDSNException
from obspy import UTCDateTime
import os

starttime = UTCDateTime("2012-10-01T01:00:00")
endtime = UTCDateTime("2012-10-02T01:00:00")

output_file = "/g/data/ha3/US_test.h5"


temp_sta = "249A"
perm_sta = "255A"


client = Client("IRIS")

ref_inv = client.get_stations(network="TA",
                              station=perm_sta,
                              channel="BHZ",
                              starttime=starttime,
                              endtime=endtime,
                              level='channel')

print(ref_inv)

temp_inv = client.get_stations(network="TA",
                               station=temp_sta,
                               channel="BHZ",
                               starttime=starttime,
                               endtime=endtime,
                               level='channel')

print(temp_inv)

# get waveforms
temp_st = client.get_waveforms(network="TA",
                               station=temp_sta,
                               channel="BHZ",
                               location="*",
                               starttime=starttime,
                               endtime=endtime)

ref_st = client.get_waveforms(network="TA",
                              station=perm_sta,
                              channel="BHZ",
                              location="*",
                              starttime=starttime,
                              endtime=endtime)

print(temp_st)

ref_st[0].stats.network = "XX"
ref_inv[0].code = "XX"

print(ref_st)
print(ref_inv)


if os.path.exists(output_file):
    os.remove(output_file)


# open up asdf
ds = pyasdf.ASDFDataSet(output_file)

ds.add_stationxml(ref_inv)
ds.add_stationxml(temp_inv)


# uniquely tag the reference station data to the temporary station data so that they can be matched together later
ds.add_waveforms(ref_st, tag="id001", labels=["region_1"])
ds.add_waveforms(temp_st, tag="raw_recording", labels=["region_1", "id001"])

print(ds.waveforms.list())
for sta_acc in ds.waveforms.list():
    print(ds.waveforms[sta_acc].list())

del ds