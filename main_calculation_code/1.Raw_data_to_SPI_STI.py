# Calculate SPI (STI) based on precipitation (air temperature) data under different scenarios

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool
from scipy import stats
import os
from scipy.stats import genextreme as gev
from scipy.stats import gamma
from scipy.stats import lognorm
from scipy.stats import norm
from scipy.stats import gennorm
from scipy.stats import expon
from scipy.stats import pearson3
import scipy.io as scio
import scipy.stats as stats
from scipy.stats import gennorm




def drought_index_jizhun(c1,c2):
    c=c1
    output=np.zeros((len(c)))
    if sum(c)>0.000000001 and sum(c==0)<=0.5*len(c)  and np.std(c)>0.000000001 and sum(c2)>0.000000001 and np.std(c2)>0.000000001 :
        shuzu_nozero=c2[c2!=0]
        xuanzehoudefenbu=selectdistribution(shuzu_nozero)
        if xuanzehoudefenbu !=1:
            t=xuanzehoudefenbu.fit(shuzu_nozero)
            if len(t)==2:
                feiling_changdu=len(shuzu_nozero)
                for i in range(0,len(c)):
                    if c[i]!=0:
                        x=c[i]
                        output[i]=xuanzehoudefenbu.cdf(x,t[0],t[1]) * feiling_changdu / (len(c)) + 1 - (feiling_changdu / (len(c)))
                    else:
                        x=c[i]
                        output[i]= 1 - (feiling_changdu / (len(c)))
            elif len(t)==3:
                feiling_changdu=len(shuzu_nozero)
                for i in range(0,len(c)):
                    if c[i]!=0:
                        x=c[i]
                        output[i]=xuanzehoudefenbu.cdf(x,t[0],t[1],t[2]) * feiling_changdu / (len(c)) + 1 - (feiling_changdu / (len(c)))
                    else:
                        x=c[i]
                        output[i]= 1 - (feiling_changdu / (len(c)))
        elif xuanzehoudefenbu==1:
            xiaoda= np.sort(c)
            xiaoda = xiaoda.tolist()
            for l in range(0,len(c)):
                output[l]=(xiaoda.index(c[l])+1)/len(c)
        for i in range(0, len(c)):
            if output[i] > 0 and output[i] <= 0.5:
                tk = (np.log(1 / (output[i] ** 2))) ** 0.5
                SPI = -(tk - (2.515517 + 0.802583 * tk + 0.010328 * (tk ** 2)) / (
                        1 + 1.432788 * tk + 0.189269 * (tk ** 2) + 0.001308 * (tk ** 3)))
                output[i] = SPI
            elif output[i] > 0.5 and output[i] < 1:
                tk = (np.log(1 / (1 - output[i]) ** 2)) ** 0.5
                SPI = (tk - (2.515517 + 0.802583 * tk + 0.010328 * (tk ** 2)) / (
                        1 + 1.432788 * tk + 0.189269 * (tk ** 2) + 0.001308 * (tk ** 3)))
                output[i] = SPI
    return output

def selectdistribution(c):
    distribution=norm,lognorm,pearson3,gamma,gev,gennorm
    for fenbu in distribution:
        t=fenbu.fit(c)
        b1 = stats.kstest(c, fenbu.cdf, args=(t), alternative='two-sided')
        if b1[1]>0.05:
            return fenbu
            break
        else:
            continue
    return 1

def fenyuejisuanganhanzhishu_jizhun(c):
    c1=c[0]
    c2=c[1]
    output=np.zeros((int(len(c[0])/12),12))
    c1=c1.reshape((int(len(c[0])/12),12))
    c2=c2.reshape((int(len(c[0])/12),12))
    for i in range(0,12):
        output[:,i]=drought_index_jizhun(c1[:,i],c2[:,i])
    shuchu=output.flatten()
    return shuchu


def bijiaojida_xiuzheng(c,d):
    output=np.zeros((len(c),len(c[0])))
    chushi=c
    c=c.reshape((len(c),int(len(c[0])/12),12))
    d=d.reshape((len(d),int(len(d[0])/12),12))
    for i in range(0,len(c)):
        for j in range(0,len(c[0][0])):
            if sum(c[i,:,j])!=0 and sum(d[i,:,j])!=0:
                dangyue_d=d[i,:,j]
                dangyue_d=dangyue_d[dangyue_d!=0]
                d_max=max(dangyue_d)
                d_min=min(dangyue_d)
                for k in range(0,len(c[0])):
                    if c[i,k,j]!=0 and c[i,k,j]>d_max:
                        output[i,12*k+j]=d_max
                    elif c[i,k,j]!=0 and c[i,k,j]<d_min:
                        output[i,12*k+j]=d_min
                    else:
                        output[i,12*k+j]=(c[i,k,j])
    return output




if __name__=="__main__":
    shijian = 'E:\\his\\'
    danbianliang='pr\\'
    his_path=r'E:\\his\\weidu_zhuanhuahou\\pr\\'
    path_list = os.listdir(his_path)
    his_path_jizhun=r'E:\\his\\weidu_zhuanhuahou\\'+danbianliang
    path_list_jizhun=os.listdir(his_path_jizhun)
    for name_id in range(0,len(path_list)):
        file_path=his_path+'\\'+path_list[name_id]
        mrro_data=np.load(file=file_path)
        mrro_data=mrro_data[:,:]
        jizhun_data=np.load(his_path_jizhun+'\\'+path_list_jizhun[name_id])
        jizhun_data=jizhun_data
        if shijian == 'E:\\his\\':
            mrro_data = mrro_data[:, 1200:1980]
            jizhun_data = jizhun_data[:, 1200:1980]
        if shijian == 'E:\\his_nat\\':
            mrro_data = mrro_data[:, 1200:1980]
            jizhun_data = jizhun_data[:, 1200:1980]
        mrro_data=np.nan_to_num(mrro_data)
        jizhun_data=np.nan_to_num(jizhun_data)
        mrro_data_yuechidu=mrro_data
        jizhun_data_yuechidu=jizhun_data
        mrro_data_yuechidu = bijiaojida_xiuzheng(mrro_data_yuechidu, jizhun_data_yuechidu)
        p = Pool(14)
        resu = []
        result = []
        for j in range(0, len(mrro_data_yuechidu)):
            c = [mrro_data_yuechidu[j, :], jizhun_data_yuechidu[j, :]]
            resu.append(p.apply_async(fenyuejisuanganhanzhishu_jizhun, args=(c,)))
        p.close()
        p.join()
        for i in resu:
            result.append(i.get())
        data_spi_non = np.array(result)
        data_spi = np.nan_to_num(data_spi_non)
        np.save(file="D:/LS3MIP/output_zhishu/his_LS/bijiao_30_newdata/danbianliang/SPI_python/{}".format(path_list[name_id]),arr=data_spi)
