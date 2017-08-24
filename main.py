from PyQt4 import QtCore, QtGui, QtWebKit, QtNetwork, uic
import pyasdf
import os
from obspy import Stream, Trace
from obspy import UTCDateTime
import numpy as np

import sys
sys.path.insert(0, "/g/data/ha3/axc547/PycharmProjects/passive-seismic")

from xcorqc.xcorqc import IntervalStackXCorr, saveXCorrPlot
from DateAxisItem import DateAxisItem


# load in Qt Designer UI files
qc_xcor_gui_ui = "qc_xcor_gui.ui"

Ui_MainWindow, QtBaseClass = uic.loadUiType(qc_xcor_gui_ui)

class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    """
        Main Window for QC XCOR GUI application
    """

    def __init__(self):
        super(MainWindow, self).__init__()
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        self.show()
        self.raise_()


    def closeEvent(self, QCloseEvent):
        # necessary to ensure data is written into ASDF file
        try:
            # close all datasets
            del self.new_ds
            del self.ds
        except AttributeError:
            # there is no loaded in ASDF data
            pass

    def on_open_pushButton_released(self):
        self.ds_filename = str(QtGui.QFileDialog.getOpenFileName(
            parent=self, caption="Choose ASDF file",
            directory=os.path.expanduser("~"),
            filter="ASDF Files (*.h5)"))
        if not self.ds_filename:
            return

        # self.ds_filename = '/g/data/ha3/US_test.h5'
        self.out_filename = os.path.join(os.path.dirname(self.ds_filename), "output_file.h5")

        if os.path.exists(self.out_filename):
            os.remove(self.out_filename)

        self.temp_net = "TA"


        # open up the input asdf file
        self.ds = pyasdf.ASDFDataSet(self.ds_filename)

        print("ASDF File Opened")

    def on_xcor_pushButton_released(self):
        print("Performing XCOR!!")
        # dictionary with unique ids as keys and stream objects (for permanent station data) as values
        id_st_dict = {}

        # first extract obspy stream data from the ASDF file that is from a permanent station (i.e. not the temporary network)
        for net_sta in self.ds.waveforms.list():
            iter_net = net_sta.split('.')[0]

            if not self.temp_net == iter_net:
                # get the station accessor
                iter_sta_accessor = self.ds.waveforms[net_sta]

                # get a list of the waveforms
                iter_waveforms = iter_sta_accessor.list()

                for wave in iter_waveforms:
                    if not wave == "StationXML":
                        iter_st = iter_sta_accessor[wave]
                        print(iter_st)

                        # add it to the dict
                        tag = wave.split('__')[3]

                        id_st_dict[tag] = iter_st

        print(id_st_dict)

        print('')

        def test_process(st, inv):

            xcor_st = Stream()

            for tr in st:
                temp_st = Stream(traces=[tr])
                print('')

                # get the uid label
                uid_label = tr.stats.asdf.labels[1]

                print(uid_label)

                perm_st = id_st_dict[uid_label]

                print(temp_st)
                print(perm_st)
                print("DO XCOR......")

                y, x, comp = IntervalStackXCorr(perm_st, temp_st, flo=0.9, fhi=1.1, interval_seconds=86400,
                                                window_seconds=1800)

                # saveXCorrPlot(y, x, '/g/data/ha3/', 'test_plot', comp)



                for day_stack_xcor in y:
                    print(type(day_stack_xcor))
                    print(day_stack_xcor[0])

                    # fill in headers
                    stats = {'network': tr.stats.network, 'station': tr.stats.station, 'location': "",
                             'channel': uid_label[2:], 'npts': len(day_stack_xcor[0]),
                             'sampling_rate': tr.stats.sampling_rate,
                             'mseed': {'dataquality': 'D'}}

                    stats["starttime"] = tr.stats.starttime

                    xcor_st += Trace(data=day_stack_xcor[0], header=stats)


                    # break

            return xcor_st

        self.ds.process(process_function=test_process,
                   output_filename=self.out_filename,
                   tag_map={"raw_recording": "xcor"})

        # load in new ASDF that has been written
        self.new_ds = pyasdf.ASDFDataSet(self.out_filename)

        # now just copy over waveforms from input ASDF file to new ASDF
        for sta_id in self.ds.waveforms.list():
            sta_acc = self.ds.waveforms[sta_id]
            for asdf_st_id in sta_acc.list():
                if not asdf_st_id == "StationXML":
                    self.new_ds.add_waveforms(sta_acc[asdf_st_id], tag=asdf_st_id.split('__')[3])
                else:
                    self.new_ds.add_stationxml(sta_acc[asdf_st_id])



        del self.ds

        self.id_st_dict = id_st_dict

        # now plot the data
        self.update_plots()

    def update_plots(self):
        print("updating plots")

        self.waveform_graph.clear()
        self.xcor_graph.clear()

        self._state = {}
        self._state["waveform_plots"] = []

        for sta_id in self.new_ds.waveforms.list():

            sta = self.new_ds.waveforms[sta_id]

            print(sta.list())

            for wave_id in sta.list():
                if wave_id == "StationXML":
                    continue

                if not "raw_recording" in wave_id:
                    continue

                raw_st = sta[wave_id]

                # get the id for the xcor and perm station
                uid_label = raw_st[0].stats.asdf.labels[1]

                # get asscoaited perm stream
                perm_st = self.id_st_dict[uid_label]

                # get the xcor
                for wave_id_2 in sta.list():
                    if wave_id_2 == "StationXML":
                        continue

                    loc_id = wave_id_2.split('__')[0].split('.')[3]
                    print(loc_id)
                    print(uid_label)
                    xcor_id = "id" + loc_id

                    if xcor_id == uid_label:
                        xcor_st = sta[wave_id_2]
                        break

                print(raw_st[0])
                print(raw_st[0].times())
                print(raw_st[0].stats.starttime.timestamp)

                waveform_plot = self.waveform_graph.addPlot(
                    0, 0, title=raw_st[0].id,
                    axisItems={'bottom': DateAxisItem(orientation='bottom',
                                                      utcOffset=0)})

                waveform_plot.plot(raw_st[0].times() + raw_st[0].stats.starttime.timestamp, raw_st[0].data, pen="b")

                self._state["waveform_plots"].append(waveform_plot)

                waveform_plot = self.waveform_graph.addPlot(
                    1, 0, title=perm_st[0].id,
                    axisItems={'bottom': DateAxisItem(orientation='bottom',
                                                      utcOffset=0)})

                waveform_plot.plot(perm_st[0].times() + perm_st[0].stats.starttime.timestamp, perm_st[0].data, pen="r")

                self._state["waveform_plots"].append(waveform_plot)

                waveform_plot.show()

                xcor_x_axis = np.linspace(-3600, 3600, xcor_st[0].stats.npts)
                xcor_plot = self.xcor_graph.addPlot(0,0)
                xcor_plot.plot(xcor_x_axis, xcor_st[0].data)

                xcor_plot.show()

        for plot in self._state["waveform_plots"][1:]:
            plot.setXLink(self._state["waveform_plots"][0])
            plot.setYLink(self._state["waveform_plots"][0])



            #     raw_wave_x_axis = np.linspace(0, raw_st[0].stats.npts, raw_st[0].stats.npts)
            #     perm_wave_x_axis = np.linspace(0, raw_st[0].stats.npts, perm_st[0].stats.npts)
            #     xcor_x_axis = np.linspace(-3600, 3600, raw_st[0].stats.npts)
            #
            #     plt.subplot(2, 2, sub_count)
            #     plt.plot(perm_wave_x_axis, perm_st[0].data, c="b")
            #
            #     plt.subplot(2, 2, sub_count)
            #
            #     plt.plot(raw_wave_x_axis, raw_st[0].data, c="r")
            #     sub_count += 1
            #
            #     plt.subplot(2, 2, sub_count)
            #     plt.plot(xcor_st[0].data)
            #
            #
            # plt.show()


if __name__ == "__main__":
    app = QtGui.QApplication([])
    w = MainWindow()
    w.raise_()
    app.exec_()