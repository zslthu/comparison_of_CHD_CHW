# Second, the new joint probability was calculated by combining the modified SPI and STI sets 
# and then compared with the original joint probability to obtain the attribution results.

import numpy as np
import pyvinecopulib as pv
import os
import re
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool
from scipy.stats import linregress

def empirical_distribution_function(values):
    sorted_indices = np.argsort(-values)
    sorted_order = np.argsort(sorted_indices) + 1
    result = (len(values) - sorted_order + 0.5) / len(values)
    return result

def raw_data_standardized_index_joint_probability(c1):
    SPI=c1[0]
    STI=c1[1]
    SPI_jizhun=c1[2]
    STI_jizhun=c1[3]
    if sum(SPI)!=0 and sum(STI)!=0 and sum(SPI_jizhun)!=0 and sum(STI_jizhun)!=0:
        joint_probability_hot_dry=calculation_joint_probability(SPI,STI,SPI_jizhun,STI_jizhun)
    else:
        joint_probability_hot_dry=np.zeros((len(SPI)))

    return joint_probability_hot_dry

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

def calculation_joint_probability(SPI,STI,SPI_jizhun,STI_jizhun):
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

    return copula_cdf_jizhun


def calculate_trends(data):
    num_grid_points, num_months = data.shape
    trends = np.zeros(num_grid_points)
    time = np.arange(num_months)
    for i in range(num_grid_points):
        slope, _, _, _, _ = linregress(time, data[i, :])
        trends[i] = slope

    return trends

if __name__=="__main__":
    SPEI_data_path='E:/his_nat/output_zhishu/his_hisnat/bijiao_30_newdata_jiyuhis/danbianliang_jiyuhis_fenyue/SPI_python_xiuzheng/'
    STI_data_path='E:/his_nat/output_zhishu/his_hisnat/bijiao_30/danbianliang/STI_python_jizhun_xiuzheng/'
    SPEI_data_jizhun_path='E:/his/output_zhishu/his_hisnat/bijiao_30/danbianliang/SPI_python_jizhun_xiuzheng/'
    STI_data_jizhun_path='E:/his/output_zhishu/his_hisnat/bijiao_30/danbianliang/STI_python_jizhun_xiuzheng/'
    model_list=['ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','EC-Earth3','FGOALS-g3','MRI-ESM2-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-LR']
    for xishu in [-1,1]:
        for model_i in range(0,len(model_list)):
            SPEI_data=np.load(SPEI_data_path+'pr'+model_list[model_i]+'.npy')
            STI_data=np.load(STI_data_path+'tas_Amon_'+model_list[model_i]+'.npy')
            SPEI_data_jizhun=np.load(SPEI_data_jizhun_path+'pr'+model_list[model_i]+'.npy')
            STI_data_jizhun=np.load(STI_data_jizhun_path+'tas_Amon_'+model_list[model_i]+'.npy')
            SPEI_data[SPEI_data>100]=0
            STI_data[STI_data>100]=0
            SPEI_data_jizhun[SPEI_data_jizhun>100]=0
            STI_data_jizhun[STI_data_jizhun>100]=0
            SPEI_data[SPEI_data<-100]=0
            STI_data[STI_data<-100]=0
            SPEI_data_jizhun[SPEI_data_jizhun<-100]=0
            STI_data_jizhun[STI_data_jizhun<-100]=0
            SPEI_data=SPEI_data*xishu
            SPEI_data_jizhun=SPEI_data_jizhun*xishu
            STI_data=-STI_data
            STI_data_jizhun=-STI_data_jizhun
            p = Pool(14)
            resu = []
            result = []
            for j in range(0, len(SPEI_data)):
                c = [SPEI_data[j,:],STI_data[j,:],SPEI_data_jizhun[j,:],STI_data_jizhun[j,:]]
                resu.append(p.apply_async(raw_data_standardized_index_joint_probability, args=(c,)))
            p.close()
            p.join()
            for i in resu:
                result.append(i.get())
            data_spi_non = np.array(result)
            data_spi = np.nan_to_num(data_spi_non)
            print(len(data_spi[data_spi < -2]))
            if xishu==1:
                np.save('E:/his_nat/output_zhishu/his_hisnat/bijiao_30_newdata_jiyuhis/lianhegailv_0903/fuhe_jiyuhis_spi_fenyue_xiuzheng/SP-TI_R/'+model_list[model_i]+'.npy',arr=data_spi)
            if xishu==-1:
                np.save('E:/his_nat/output_zhishu/his_hisnat/bijiao_30_newdata_jiyuhis/lianhegailv_0903/fuhe_jiyuhis_spi_fenyue_xiuzheng/-SP-TI_R/'+model_list[model_i]+'.npy',arr=data_spi)
