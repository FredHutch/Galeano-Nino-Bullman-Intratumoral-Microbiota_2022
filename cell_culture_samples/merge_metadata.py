#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import sys
pd.set_option('precision', 0)

def merge_cellsmeta2(df1,df2):
    df_merged = pd.concat([df1, df2], sort=False)
    df_merged = df_merged.round()
    return df_merged

def feed_csvs(path):#this will return a list of csvs in your path
    file_list = os.listdir(path)
    csv_list = []
    for each_file in file_list:
        if each_file.endswith('genus.csv'):
        #if each_file.endswith('merged.csv'):
            csv_list.append(path+'/'+each_file)
    return csv_list

if __name__ == "__main__":
    path = sys.argv[1]
    csv_merged = path+'/csv_novami.csv'
    csv_list = feed_csvs(path)

    csv1 = csv_list[0]
    df1 = pd.read_csv(csv1,header = 0,sep = ',')

    for each_csv in csv_list[1:]:
        print(each_csv)
        df2 = pd.read_csv(each_csv,header = 0,sep = ',')
        df1 = merge_cellsmeta2(df1,df2)

    #print(df1)
    df1 = df1.fillna(0)
    df1.to_csv(csv_merged,sep=',',index=False)

