#!/usr/bin/env python
import numpy as np
import pandas as pd

def main():
    # # Allowing user to specify the import everytime
    data_path = input("Enter file name of List of samples in order they were put into QPCR: ")
    print(data_path)
    data_list = pd.read_csv(f'{data_path}', header=None)
    data_list.to_numpy
    # # making new array that is 3x the length
    data_listX3 = np.zeros((len(data_list)*3, 1),  dtype=object)

    # # Filling in new array
    for i in range(len(data_list)):
        top_num = (i+1)*3 - 1
        mid_num = top_num - 1
        bot_num = mid_num - 1
        data_listX3[top_num][0] = data_list[0][i]
        data_listX3[mid_num][0] = data_list[0][i]
        data_listX3[bot_num][0] = data_list[0][i]

    np.savetxt("data_listX3.csv", data_listX3, fmt='%s')
main()