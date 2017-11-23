class DustIndex:
    from netCDF4 import Dataset
    import pandas as pd
    import numpy as np
    from datetime import datetime as dt
    from datetime import timedelta

    def __init__(self,datestr,NC_PATH):
        self.dataFrame = self.__openNcFile(datestr,NC_PATH)
        self.DI_low,self.DI_high,self.updateTime = self.__getDI()

        if self.DI_high > 0.6:
            self.DI_high = 0.6

        if self.DI_low > 0.6:
            self.DI_low = 0.6


    def __openNcFile(self,datestr,NC_PATH):
        # Parameter:
        NC_NAME = "/li{}.b532".format(datestr[2:])

        # Get data from nc-file:
        print(datestr)
        print(NC_PATH+NC_NAME)
        nc = self.Dataset(NC_PATH + NC_NAME)
        dustIndexLow = nc.variables["DustIndexLowLayer"][:].copy()
        dustIndexTotal = nc.variables["DustIndexTotal"][:].copy()
        seconds = nc.variables["Time"][:].copy()
        nc.close()

        # Handle missing Values:
        dustIndexTotal[self.np.where(dustIndexTotal > 1e30)] = self.np.nan
        dustIndexLow[self.np.where(dustIndexLow > 1e30)] = self.np.nan

        dustIndexHigh = self.np.subtract(dustIndexTotal,dustIndexLow)

        # Shape time to the right format:
        time = []
        start_time = self.dt(int(datestr[:4]), int(datestr[4:6]), int(datestr[6:]), 0, 0, 0)
        for t in seconds:
            time.append(start_time + self.timedelta(seconds=int(t)))
        time = self.np.asarray(time)

        df = self.pd.DataFrame({'time': time, 'DIL': dustIndexLow, 'DIH': dustIndexHigh})
        # df = df[::-1]  # reverse columns
        return df

    def __getDI(self):
        __dil = self.dataFrame["DIL"][-30:].mean() # Each two minutes we get a measurement. --> last hour = the last 30 values
        __dih = self.dataFrame["DIH"][-30:].mean()
        __up_time = self.dataFrame["time"].tail(1).values[0]
        __ts = self.pd.to_datetime(str(__up_time))
        __d = __ts.strftime('%Y %m %d, %H:%M:%S')
        return(__dil,__dih,__d)