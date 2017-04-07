# phys_toolkit.py
#
# David Lister
# 2017-02-18
#

import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.odr as odr
from uncertainties import ufloat
import math

class CsvInfo:
    def __init__(self, fname, x, y, x_conv_factor, y_conv_factor):
        self.fname = fname
        self.x = x
        self.y = y
        self.x_conv_factor = x_conv_factor
        self.y_conv_factor = y_conv_factor

class CsvRawData:
    def __init__(self, names, data, errors):
        self.names = names
        self.data = data
        self.error = error

class Data:
    def __init__(self, x_lst, y_lst, x_err, y_err):
        self.x_lst = x_lst
        self.y_lst = y_lst
        self.x_err = x_err
        self.y_err = y_err

class DataCollection:
    def __init__(self, name, csvInfo = None, csvData = None, odrParams = None, odrFit = None, conversions = None, data = None, plotInfo = None):
        self.name = name
        self.csvInfo = csvInfo
        self.csvData = csvData
        self.odrParams = ordParams
        self.odrFit = odrFit
        self.conversions = conversions
        self.data = data

    def openCSV(self):
        self.csvData = CsvRawData(csv_to_lst_with_error(self.csv.fname))

    def processData(self):
        self.data = Data(convert_data_and_error(self.csvData.data, self.csvData.error, self.csvInfo.x, self.csvInfo.y, self.csvInfo.x_conv_factor, self.csvInfo.y_conv_factor))

    def odrFit(self):
        pass
        


def csv_to_lst_with_error(fname):
    with open(fname, 'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(4096))
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)
        out = []
        for row in reader:
            out.append(row)

        cols = len(out[-1])
        error = out[0:cols]
        names = out[cols]
        out = out[cols+1:]

        error = np.array(error)
        error = error.astype(float)
        names = np.array(names)
        out = np.array(out)
        out = out.astype(float)
        out = out.transpose()
        return names, out, error

    
def linear(p, x):
    m, b = p
    return m*x + b

def sine_model(p, x):
    A, f, phi, b = p
    return A * np.sin( 2 * np.pi * f * x - phi) + b

linearModel = odr.Model(linear)
sineModel = odr.Model(sine_model)

one_to_one = lambda x: x

gauss_to_tesla = lambda x: x/10000.0



def convert_data_and_error(data, err_params, x, y, x_conv_factor, y_conv_factor):
    x_lst = []
    y_lst = []
    x_err = []
    y_err = []
    for i in range(len(data[x])):
        x_temp = ufloat(data[x][i], abs(err_params[x][0] * data[x][i]) + err_params[x][1])
        x_temp = x_conv_factor(x_temp)
        x_lst.append(x_temp.nominal_value)
        x_err.append(x_temp.std_dev)
        y_temp = ufloat(data[y][i], abs(err_params[y][0] * data[y][i]) + err_params[y][1])
        y_temp = y_conv_factor(y_temp)
        y_lst.append(y_temp.nominal_value)
        y_err.append(y_temp.std_dev)

    return x_lst, y_lst, x_err, y_err

def fit_csv_data_to_model(fname, x, y, x_conv_factor, y_conv_factor, model):
    names, data, err_params = csv_to_lst_with_error(fname)
    x_lst, y_lst, x_err, y_err = convert_data_and_error(data, err_params, x, y, x_conv_factor, y_conv_factor)

    odr_data = odr.RealData(x_lst, y_lst, x_err, y_err)
    odr_model = odr.ODR(odr_data, model, beta0=est)
    fit = odr_model.run()

    out = []
    for i in range(len(fit.beta)):
        out.append(ufloat(fit.beta[i], fit.sd_beta[i]))

    return odr_data, out, fit


def fit_data_to_model(x, y, x_err, y_err, model, est):
    odr_data = odr.RealData(x, y, x_err, y_err)
    odr_model = odr.ODR(odr_data, model, beta0=est)
    fit = odr_model.run()

    out = []
    for i in range(len(fit.beta)):
        out.append(ufloat(fit.beta[i], fit.sd_beta[i]))

    return odr_data, out, fit
