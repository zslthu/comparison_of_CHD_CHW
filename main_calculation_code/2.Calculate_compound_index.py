# Calculate compound hot-dry index and compound hot-wet index using SPI and STI under different scenarios

import numpy as np
from scipy.stats import genextreme as gev
from scipy.stats import gamma,lognorm,norm,gennorm,expon,pearson3
from scipy import stats
import pyvinecopulib as pv
import matplotlib.pyplot as plt
import os
import re
from scipy.io import loadmat
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import linregress

def empirical_distribution_function(values):
    sorted_indices = np.argsort(-values)
    sorted_order = np.argsort(sorted_indices) + 1
    result = (len(values) - sorted_order + 0.5) / len(values)
    return result

def raw_data_standardized_index_composite_index(c1):
    SPI=c1[0]
    STI=c1[1]
    SPI_jizhun=c1[2]
    STI_jizhun=c1[3]
    if sum(SPI)!=0 and sum(STI)!=0 and sum(SPI_jizhun)!=0 and sum(STI_jizhun)!=0:
        composite_index_hot_dry=calculation_composite_index(SPI,STI,SPI_jizhun,STI_jizhun)
    else:
        composite_index_hot_dry=np.zeros((len(SPI)))

    return composite_index_hot_dry

def index_paixu_jizhun(value, jizhun1):
    output = np.zeros(len(value))
    for i in range(len(value)):
        combined = np.append(value[i], jizhun1)
        rank = (len(combined) - np.sum(combined >= value[i]) + 0.5) / len(combined)
        output[i] = rank
        if rank < 0:
            print(rank)
        if rank > 1:
            print(rank)
    return output

def calculation_composite_index(SPI,STI,SPI_jizhun,STI_jizhun):
    SPI_data_em=index_paixu_jizhun(SPI,SPI_jizhun)
    STI_data_em=index_paixu_jizhun(STI,STI_jizhun)
    combined_array = np.column_stack((SPI_data_em, STI_data_em))
    copula_set = [pv.BicopFamily.gaussian, pv.BicopFamily.student, pv.BicopFamily.clayton
                , pv.BicopFamily.gumbel, pv.BicopFamily.frank, pv.BicopFamily.joe]
    min_bic = float('inf')
    best_copula = None
    for copula in copula_set:
        cop = pv.Bicop(copula)
        bic = cop.bic(combined_array)
        if bic < min_bic:
            min_bic = bic
            best_copula = copula
    copula = pv.Bicop(best_copula)
    copula.fit(combined_array)
    copula_cdf=copula.cdf(combined_array)
    SPI_data_em_jizhun=index_paixu_jizhun(SPI_jizhun,SPI_jizhun)
    STI_data_em_jizhun=index_paixu_jizhun(STI_jizhun,STI_jizhun)
    combined_array_jizhun = np.column_stack((SPI_data_em_jizhun, STI_data_em_jizhun))
    copula_set_jizhun = [pv.BicopFamily.gaussian, pv.BicopFamily.student, pv.BicopFamily.clayton
                , pv.BicopFamily.gumbel, pv.BicopFamily.frank, pv.BicopFamily.joe]
    min_bic = float('inf')
    best_copula = None
    for copula_jizhun in copula_set:
        cop = pv.Bicop(copula_jizhun)
        bic = cop.bic(combined_array_jizhun)
        if bic < min_bic:
            min_bic = bic
            best_copula_jizhun = copula_jizhun
    copula_jizhun = pv.Bicop(best_copula_jizhun)
    copula_jizhun.fit(combined_array_jizhun)
    copula_cdf_jizhun=copula_jizhun.cdf(combined_array_jizhun)
    output=index_paixu_jizhun(copula_cdf,copula_cdf_jizhun)
    for i in range(0, len(output)):
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



if __name__=="__main__":
    SPI_data_path='E:\\his_nat\\output_zhishu\\his_his_nat\\bijiao_30\\danbianliang\\STI_python_jizhun_xiuzheng_new\\'
    STI_data_path='E:\\his_nat\\output_zhishu\\his_his_nat\\bijiao_30\\danbianliang\\STI_python_jizhun_xiuzheng_new\\'
    SPI_data_jizhun_path='E:\\his\\output_zhishu\\his_his_nat\\bijiao_30\\danbianliang\\STI_python_jizhun_xiuzheng_new\\'
    STI_data_jizhun_path='E:\\his\\output_zhishu\\his_his_nat\\bijiao_30\\danbianliang\\STI_python_jizhun_xiuzheng_new\\'
    model_list=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','EC-Earth3','FGOALS-g3','MRI-ESM2-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-LR']
    for xishu in [-1,1]:
        for model_i in range(0,len(model_list)):
            SPI_data=np.load(SPI_data_path+'pr'+model_list[model_i]+'.npy')
            STI_data=np.load(STI_data_path+'tas_Amon_'+model_list[model_i]+'.npy')
            SPI_data_jizhun=np.load(SPI_data_jizhun_path+'pr'+model_list[model_i]+'.npy')
            STI_data_jizhun=np.load(STI_data_jizhun_path+'tas_Amon_'+model_list[model_i]+'.npy')
            SPI_data[SPEI_data>100]=0
            STI_data[STI_data>100]=0
            SPI_data_jizhun[SPEI_data_jizhun>100]=0
            STI_data_jizhun[STI_data_jizhun>100]=0
            SPI_data[SPEI_data<-100]=0
            STI_data[STI_data<-100]=0
            SPI_data_jizhun[SPEI_data_jizhun<-100]=0
            STI_data_jizhun[STI_data_jizhun<-100]=0
            SPI_data=SPI_data*xishu
            SPI_data_jizhun=SPI_data_jizhun*xishu
            STI_data=-STI_data
            STI_data_jizhun=-STI_data_jizhun
            p = Pool(14)
            resu = []
            result = []
            for j in range(0, len(SPEI_data)):
                c = [SPI_data[j,:],STI_data[j,:],SPI_data_jizhun[j,:],STI_data_jizhun[j,:]]
                resu.append(p.apply_async(raw_data_standardized_index_composite_index, args=(c,)))
            p.close()
            p.join()
            for i in resu:
                result.append(i.get())
            data_spi_non = np.array(result)
            data_spi = np.nan_to_num(data_spi_non)
            print(len(data_spi[data_spi < -2]))
            if xishu==1:
                np.save('D:\\chd_chw_SPI\\fuhe\\hisnat_1980_2014\\SP-TI\\SP-TI'+model_list[model_i]+'.npy',arr=data_spi)
            if xishu==-1:
                np.save('D:\\chd_chw_SPI\\fuhe\\hisnat_1980_2014\\-SP-TI\\-SP-TI'+model_list[model_i]+'.npy',arr=data_spi)
