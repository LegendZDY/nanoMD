#!/usr/bin/env python3
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20201224
Description:
This Script of ONTrnaDirectSeq from Legendzdy.  Of course, this means there's a possibility
for other ways.Use at your own discretion.
"""
import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
import joblib, gzip
import argparse


def fetch_reads(sample, model):
    Snames = ["read", "UmisRate", "UmisQrate", "t_a_rate", "t_a_Qrate",
              "t_c_rate", "t_c_Qrate", "t_g_rate", "t_g_Qrate", "del_rate"]
    Sdf = pd.read_csv(sample)
    Sdf = Sdf[Snames]
    Sdf = Sdf.replace([np.inf, -np.inf], np.nan)
    Sdf.eval("""
    Um_UmQ = UmisRate*UmisQrate
    Um_t_a_rate = UmisRate*t_a_rate
    Um_t_a_Qrate = UmisRate*t_a_Qrate
    Um_t_c_rate = UmisRate*t_c_rate
    Um_t_c_Qrate = UmisRate*t_c_Qrate
    Um_t_g_rate = UmisRate*t_g_rate
    Um_t_g_Qrate = UmisRate*t_g_Qrate
    Um_del_rate = UmisRate*del_rate
    UQ_t_a_rate = UmisQrate*t_a_rate
    UQ_t_a_Qrate = UmisQrate*t_a_Qrate
    UQ_t_c_rate = UmisQrate*t_c_rate
    UQ_t_c_Qrate = UmisQrate*t_c_Qrate
    UQ_t_g_rate = UmisQrate*t_g_rate
    UQ_t_g_Qrate = UmisQrate*t_g_Qrate
    UQ_del_rate = UmisQrate*del_rate
    t_a_rate_t_a_Qrate = t_a_rate*t_a_Qrate
    t_a_rate_t_c_rate = t_a_rate*t_c_rate
    t_a_rate_t_c_Qrate = t_a_rate*t_c_Qrate
    t_a_rate_t_g_rate = t_a_rate*t_g_rate
    t_a_rate_t_g_Qrate = t_a_rate*t_g_Qrate
    t_a_rate_del_rate = t_a_rate*del_rate
    t_a_Qrate_t_c_rate = t_a_Qrate*t_c_rate
    t_a_Qrate_t_c_Qrate = t_a_Qrate*t_c_Qrate
    t_a_Qrate_t_g_rate = t_a_Qrate*t_g_rate
    t_a_Qrate_t_g_Qrate = t_a_Qrate*t_g_Qrate
    t_a_Qrate_del_rate = t_a_Qrate*del_rate
    t_c_rate_t_c_Qrate = t_c_rate*t_c_Qrate
    t_c_rate_t_g_rate = t_c_rate*t_g_rate
    t_c_rate_t_g_Qrate = t_c_rate*t_g_Qrate
    t_c_rate_del_rate = t_c_rate*del_rate
    t_c_Qrate_t_g_rate = t_c_Qrate*t_g_rate
    t_c_Qrate_t_g_Qrate = t_c_Qrate*t_g_Qrate
    t_c_Qrate_del_rate = t_c_Qrate*del_rate
    t_g_rate_t_g_Qrate = t_g_rate*t_g_Qrate
    t_g_rate_del_rate = t_g_rate*del_rate
    t_g_Qrate_del_rate = t_g_Qrate*del_rate
    """, inplace=True)
    S_test = Sdf.iloc[:, 1:]
    imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
    S_test_t = imputer.fit_transform(S_test)
    clf2 = joblib.load(model)
    S_y = clf2.predict(S_test_t)
    dict, i = {}, 0
    while i < len(S_y):
        if int(S_y[i]) == 1:
            dict[Sdf["read"][i]] = 1
        i += 1
    return dict

def new_fq(sample, model, rawfq, newfq):
    dict_new = fetch_reads(sample, model)
    with open(newfq, 'w') as fo:
        with gzip.GzipFile(rawfq, "rb") as fi:
            while True:
                try:
                    line1 = next(fi).decode().strip()
                    line2 = next(fi).decode().strip()
                    line3 = next(fi).decode().strip()
                    line4 = next(fi).decode().strip()
                    if dict_new.get(line1.split()[0][1:], None) != None:
                        fo.write("{}\n{}\n{}\n{}\n".format(
                            line1, line2, line3, line4))
                except StopIteration:
                    break

def main():
    parser = argparse.ArgumentParser(
        description='Extract length and quality information of sequences form fastq file')
    parser.add_argument('-s', '--signal', type=str, required=False,
                        help='sample raw data fecth raw signal')
    parser.add_argument('-r', '--rawfq', type=str, required=True,
                        help='sample raw fastq data')
    parser.add_argument('-m', '--model', type=str, required=True,
                        help='train model')
    parser.add_argument('-o', '--out_file', type=str, required=True,
                        help='Path to an output file to be created,default= ./')
    args = parser.parse_args()
    new_fq(args.signal, args.model, args.rawfq, args.out_file)


if __name__ == '__main__':
    main()
# new_fq("OV-FTO-2-0603.fastq.sam.csv", "test.pkl", "all.fastq", "new.fastq")
