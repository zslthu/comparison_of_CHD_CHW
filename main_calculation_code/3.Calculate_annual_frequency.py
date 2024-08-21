# Calculate the annual frequency of compound extreme events under different scenarios 

import numpy as np
import os

def andianhuatu(c,fenweishu):
    m = 0
    n = 0
    year = 0
    nianfen=int(len(c[0])/12)
    junzhi_anshijian=np.zeros((nianfen))
    junzhi_anshijian_lishi = np.zeros((nianfen))
    q10_c = np.percentile(c, fenweishu)
    while year < nianfen:
        m = n + 12
        data_year = c[:, n:m]
        tr = 0
        quanqiu = []
        quanqiu_lishi = []
        mianji_zongshumu=0
        for j in range(0,len(data_year)):
            if sum(data_year[j,:])!=0:
                nt=0
                zhibiao=[]
                f=data_year[j,:]
                for o in range(0, 12):
                    if f[o] <q10_c:
                        nt = nt + 1
                        zhibiao.append(f[o])
                if nt>=2:
                    c_year_liedu=1
                else:
                    c_year_liedu=0
                if nt!=0:
                    c_year_lishi=nt
                    tr = tr + 1
                else:
                    c_year_lishi = 0
                    tr = tr + 0
                quanqiu.append(c_year_liedu)
                quanqiu_lishi.append(c_year_lishi)
                mianji_zongshumu=mianji_zongshumu+1
            else:
                tr=tr+0
        n=m
        bu = sum(np.array(quanqiu))
        bu_lishi = sum(np.array(quanqiu_lishi))
        quanqiupingjun=bu/mianji_zongshumu
        quanqiupingjun_lishi = bu_lishi / tr
        junzhi_anshijian[year]=quanqiupingjun
        junzhi_anshijian_lishi[year] = quanqiupingjun_lishi/12
        year=year+1
    return junzhi_anshijian,junzhi_anshijian_lishi

if __name__=="__main__":
    shiqi = ['his']
    zhishu = ['-SP-TI', 'SP-TI']
    for shijian in shiqi:
        for duiyingzhishu in zhishu:
            path = 'E:\\' + shijian + '\\output_zhishu\\his_LS\\bijiao_30\\fuhe_paixu_fan\\' + duiyingzhishu+'_jizhun_xiuzheng\\'
            bapchun_path = 'E:\\' + shijian + '\\output_zhishu\\\\his_LS\\bijiao_30\\tu_paixu_fan\\zongtu\\ganhantezheng\\' + duiyingzhishu+'\\'
            path_list=os.listdir(path)
            for name in path_list:
                data=np.load(path+name)
                data=data
                mianji,pinlv=andianhuatu(data,10)
                print(mianji.shape)
                np.save(file=bapchun_path+'\\mianji\\'+name,arr=mianji)
                np.save(file=bapchun_path+'\\pinlv\\'+name,arr=pinlv)