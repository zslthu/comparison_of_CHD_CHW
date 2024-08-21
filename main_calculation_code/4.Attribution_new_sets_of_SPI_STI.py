# Analyzing the impacts of (i) precipitation change, (ii) temperature change, and (iii) changes in precipitation-temperature correlation on compound events 
# using the method from Bevacqua et al. (2019, DOI: 10.1126/sciadv.aaw5531).
# First, three new sets of SPI and STI were constructed, each modifying only one of the three factors mentioned above.

import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool
from scipy import stats
import os
from scipy.stats import genextreme as gev
from scipy.stats import gamma
from scipy.stats import lognorm
from scipy.stats import norm
from scipy.stats import gennorm
from scipy.stats import pearson3


def drought_index_jizhun(c1,c2):
    c=c1
    output=np.zeros((len(c)))
    if sum(c)!=0 and sum(c==0)<=0.5*len(c)  and np.std(c)!=0 and sum(c2)!=0 and np.std(c2)!=0 :
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

def empirical_cdf(data):
    sorted_indices = np.argsort(data)
    sorted_data = data[sorted_indices]
    cumulative_prob = np.arange(1, len(data) + 1) / len(data)
    output = np.zeros(len(data))
    output[sorted_indices] = cumulative_prob
    return output

def find_matching_positions(a, b):
    output=np.zeros(len(a))
    for i in range(0,len(output)):
        output[i] = np.argwhere(b == a[i])
    return output


def transform_data(a, b):
    a_cdf = empirical_cdf(a)
    b_cdf = empirical_cdf(b)
    shunxun=find_matching_positions(a_cdf,b_cdf)
    new_a=np.zeros(len(a))
    for i in range(0,len(a)):
        weizhi=shunxun[i]
        new_a[i]=b[int(weizhi)]
    return new_a


def quanqiujingyanleijigailvzhuanhua(a,a_jizhun):
    output=np.zeros((len(a),len(a[0])))
    for i in range(0,len(a)):
        if sum(a[i,:])!=0 and sum(a_jizhun[i,:])!=0:
            output[i,:]=fenyuetihuanleijigailv_jizhun(a[i,:],a_jizhun[i,:])
    return output

def fenyuetihuanleijigailv_jizhun(c1,c2):
    changdu=len(c1)
    output=np.zeros((int(changdu/12),12))
    c1=c1.reshape((int(changdu/12),12))
    c2=c2.reshape((int(changdu/12),12))
    for i in range(0,12):
        output[:,i]=transform_data(c1[:,i],c2[:,i])
    shuchu=output.flatten()
    return shuchu


def bijiaojida_xiuzheng(c,d):
    output=np.zeros((len(c),len(c[0])))
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
    bianliang_zong=['tas','pr']
    shuchu_zong=['STI_python_jizhun_xiuzheng_new','SPI_python_jizhun_xiuzheng_new']
    moshi_zong=['E:\\his\\']
    for jj in range(0,2):
        bianliang=bianliang_zong[jj]
        shuchu=shuchu_zong[jj]
        for shijian in moshi_zong:
            if shijian=='E:\\his\\':
                his_path=r'E:\\his\\duibiLS\\weidu_zhuanhuahou\\'+bianliang
            elif shijian=='D:\\LS3MIP\\':
                his_path=r'D:\\LS3MIP\\weidu_zhuanhuahou\\'+bianliang
            elif shijian=='E:\\his_nat\\':
                his_path=r'E:\\his_nat\\weidu_zhuanhuahou\\'+bianliang
            path_list = os.listdir(his_path)
            if bianliang=='tas':
                path_list=['tas_Amon_ACCESS-CM2.npy', 'tas_Amon_ACCESS-ESM1-5.npy', 'tas_Amon_BCC-CSM2-MR.npy',
                           'tas_Amon_CanESM5.npy', 'tas_Amon_CESM2.npy', 'tas_Amon_CNRM-CM6-1.npy',
                           'tas_Amon_E3SM-2-0.npy', 'tas_Amon_FGOALS-g3.npy', 'tas_Amon_GFDL-CM4.npy',
                           'tas_Amon_GFDL-ESM4.npy', 'tas_Amon_IPSL-CM6A-LR.npy', 'tas_Amon_MIROC6.npy',
                           'tas_Amon_MRI-ESM2-0.npy', 'tas_Amon_NorESM2-LM.npy']
            elif bianliang=='pr':
                path_list=['pr_Amon_ACCESS-CM2.npy', 'pr_Amon_ACCESS-ESM1-5.npy', 'pr_Amon_BCC-CSM2-MR.npy',
                           'pr_Amon_CanESM5.npy', 'pr_Amon_CESM2.npy', 'pr_Amon_CNRM-CM6-1.npy',
                           'pr_Amon_E3SM-2-0.npy', 'pr_Amon_FGOALS-g3.npy', 'pr_Amon_GFDL-CM4.npy',
                           'pr_Amon_GFDL-ESM4.npy', 'pr_Amon_IPSL-CM6A-LR.npy', 'pr_Amon_MIROC6.npy',
                           'pr_Amon_MRI-ESM2-0.npy', 'pr_Amon_NorESM2-LM.npy']
            if shijian=='D:\\LS3MIP\\' and bianliang=='tas':
                path_list=['tas_Amon_CESM2_amip-lfmip-pdLC_r1i1p1f1_gn_197001-210012.npy',
                           'tas_Amon_CMCC-ESM2_amip-lfmip-pdLC_r1i1p1f1_gn_198001-210012.npy',
                           'tas_Amon_CNRM-CM6-1_amip-lfmip-pdLC_r1i1p1f2_gr_198001-201412.npy',
                           'tas_Amon_EC-Earth3_amip-lfmip-pdLC_r1i1p1f1_gn_198001-210012.npy',
                           'tas_Amon_IPSL-CM6A-LR_amip-lfmip-pdLC_r1i1p1f1_gr_198001-210012.npy',
                           'tas_Amon_MIROC6_amip-lfmip-pdLC_r1i1p1f1_gn_198001-210012.npy',
                           'tas_Amon_MPI-ESM1-2-LR_amip-lfmip-pdLC_r1i1p1f1_gn_198001-210012.npy']
            elif shijian=='D:\\LS3MIP\\' and bianliang=='pr':
                path_list=['pr_Amon_CESM2_amip-lfmip-pdLC_r1i1p1f1_gn_197001-210012.npy',
                           'pr_Amon_CMCC-ESM2_amip-lfmip-pdLC_r1i1p1f1_gn_198001-210012.npy',
                           'pr_Amon_CNRM-CM6-1_amip-lfmip-pdLC_r1i1p1f2_gr_198001-201412.npy',
                           'pr_Amon_EC-Earth3_amip-lfmip-pdLC_r1i1p1f1_gn_198001-210012.npy',
                           'pr_Amon_IPSL-CM6A-LR_amip-lfmip-pdLC_r1i1p1f1_gr_198001-210012.npy',
                           'pr_Amon_MIROC6_amip-lfmip-pdLC_r1i1p1f1_gn_198001-210012.npy',
                           'pr_Amon_MPI-ESM1-2-LR_amip-lfmip-pdLC_r1i1p1f1_gn_198001-210012.npy']
            if bianliang=='pr':
                path_list1=['pr_Amon_ACCESS-CM2.npy', 'pr_Amon_ACCESS-ESM1-5.npy', 'pr_Amon_BCC-CSM2-MR.npy',
                           'pr_Amon_CanESM5.npy', 'pr_Amon_CESM2.npy', 'pr_Amon_CNRM-CM6-1.npy',
                           'pr_Amon_E3SM-2-0.npy', 'pr_Amon_FGOALS-g3.npy', 'pr_Amon_GFDL-CM4.npy',
                           'pr_Amon_GFDL-ESM4.npy', 'pr_Amon_IPSL-CM6A-LR.npy', 'pr_Amon_MIROC6.npy',
                           'pr_Amon_MRI-ESM2-0.npy', 'pr_Amon_NorESM2-LM.npy']
            elif bianliang=='tas':
                path_list1=['tas_Amon_ACCESS-CM2.npy', 'tas_Amon_ACCESS-ESM1-5.npy', 'tas_Amon_BCC-CSM2-MR.npy',
                           'tas_Amon_CanESM5.npy', 'tas_Amon_CESM2.npy', 'tas_Amon_CNRM-CM6-1.npy',
                           'tas_Amon_E3SM-2-0.npy', 'tas_Amon_FGOALS-g3.npy', 'tas_Amon_GFDL-CM4.npy',
                           'tas_Amon_GFDL-ESM4.npy', 'tas_Amon_IPSL-CM6A-LR.npy', 'tas_Amon_MIROC6.npy',
                           'tas_Amon_MRI-ESM2-0.npy', 'tas_Amon_NorESM2-LM.npy']
            for name_id in range(0,len(path_list)):
                file_path=his_path+'\\'+path_list[name_id]
                mrro_data=np.load(file=file_path)
                jizhun_data=np.load('E:\\his_nat\\weidu_zhuanhuahou\\'+bianliang+'\\'+path_list1[name_id])
                if path_list[name_id]=='tas_Amon_CESM2_amip-lfmip-pdLC_r1i1p1f1_gn_197001-210012.npy':
                    mrro_data=mrro_data[:,120:1572]
                if path_list[name_id]=='pr_Amon_CESM2_amip-lfmip-pdLC_r1i1p1f1_gn_197001-210012.npy':
                    mrro_data=mrro_data[:,120:1572]
                if shijian=='D:\\LS3MIP\\':
                    mrro_data=mrro_data[:,0:420]
                if shijian=='E:\\his\\':
                    mrro_data=mrro_data[:,1560:1980]
                if shijian=='E:\\his_nat\\':
                    mrro_data=mrro_data[:,1560:1980]
                jizhun_data=jizhun_data[:,1560:1980]
                new_his_data=quanqiujingyanleijigailvzhuanhua(mrro_data,jizhun_data)
                np.save(file=shijian+'output_zhishu\\his_hisnat\\weidu_zhuanhuahou_jiyuhisnat_fenyue\\'+bianliang+'\\'+path_list[name_id], arr=new_his_data)
    his_path=r'E:\\his\\output_zhishu\\his_hisnat\\weidu_zhuanhuahou_jiyuhisnat_fenyue\\'+bianliang
    path_list = os.listdir(his_path)
    his_path_jizhun=r'E:\\his\\duibiLS\\weidu_zhuanhuahou\\'+bianliang
    path_list_jizhun=os.listdir(his_path_jizhun)
    for name_id in range(8,14):
        file_path=his_path+'\\'+path_list[name_id]
        mrro_data=np.load(file=file_path)
        mrro_data=mrro_data[:,:]
        jizhun_data=np.load(his_path_jizhun+'\\'+path_list_jizhun[name_id])
        if shijian == 'E:\\his_nat\\':
            jizhun_data = jizhun_data[:, 1560:1980]
            print(shijian, jizhun_data.shape)
        mrro_data_yuechidu=mrro_data
        jizhun_data_yuechidu=jizhun_data
        mrro_data_yuechidu = bijiaojida_xiuzheng(mrro_data_yuechidu, jizhun_data_yuechidu)
        p = Pool(14)
        resu=[]
        result=[]
        for j in range(0, len(mrro_data_yuechidu)):
            c = [mrro_data_yuechidu[j, :], jizhun_data_yuechidu[j, :]]
            resu.append(p.apply_async(fenyuejisuanganhanzhishu_jizhun, args=(c,)))
        p.close()
        p.join()
        for i in resu:
            result.append(i.get())
        data_spi_non = np.array(result)
        data_spi=np.nan_to_num(data_spi_non)
        np.save(file="E:/his/output_zhishu/his_hisnat/bijiao_30_newdata_jiyuhis/danbianliang_xiangguan_fenyue/+"shuchu"+/{}".format(path_list[name_id]), arr=data_spi)