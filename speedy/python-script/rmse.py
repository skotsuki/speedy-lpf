"""
printing rmse.

usage:
  python3 rmse.py


"""

#import xarray as xr  #< TBD
import netCDF4 as nc
import os
import sys
import numpy as np
import subprocess
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import matplotlib.cm as cm

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

DATA_PATH="../DATA/"
nature_dir=DATA_PATH+"nature/"
control_filename="yyyymmddhh.ctl"

RMSE_SKIP=30*4 #Skip first 1 month
level=3#rmse level

#実験ごとのクラス。　RMSEの計算や様々なマップの作成が含まれる
class exp_result:
    #system directory, title of experiment, nature
    def __init__(self,sysdir,title,nature):
        self._sysdir=sysdir
        self.title=title
        self.meanpath=DATA_PATH+"raob/"+self._sysdir+"/anal/mean/"
        self.kldvpath=DATA_PATH+"raob/"+self._sysdir+"/gues/kldv/"
        self.nature=nature        

    #get mean variables
    def get_mean(self):
        self.mean=get_nc(self.meanpath)
        return self.mean

    #get kldv variables
    def get_kldv(self):
        self.kldv=get_nc(self.kldvpath)
        return self.kldv

    #calculate kldv mean maps
    def calc_kldvmean_map(self):
        self.meankldv_map=np.mean(self.kldv.variables['t'][RMSE_SKIP:,level,:,:],axis=0)
        print(self.meankldv_map.shape)
        return self.meankldv_map

    #rmse calculation
    def calc_rmse(self):
        self.rmselist=[]
        for i,x in enumerate(self.mean.variables['time']):
            self.rmselist.append(rmse(i,level,nature,self.mean))
        return self.rmselist

    #draw rmse timeline graph.
    def drawrmse(self):
        plt.plot(self.rmselist,label="RMSE of "+self.title)

    #draw rmse maps.
    def calc_rmsemean_map(self):
        self.meanrmse_map=np.sqrt(np.mean(
            (self.mean.variables['t'][RMSE_SKIP:,level,:,:]-
             self.nature.variables['t'][RMSE_SKIP:-2,level,:,:])**2
            ,axis=0))
        return self.meanrmse_map;

    #calculate mean tempurture map
    def get_meant_map(self):
        self.meant_map=np.mean(self.mean.variables['t'][RMSE_SKIP:,level,:,:],axis=0)
        return self.meant_map

    #calculate mean of rmse map.
    def calc_meanrmse(self):
        return np.mean(self.rmselist[RMSE_SKIP:]) #skip 1 month

    # plot first temp
    def plott(self,i,ax):
        imshow_map(self.mean.variables['t'][0][level],self.title,ax)

    #plot mean temp map
    def plottmean(self,i,ax):
        imshow_map(self.meant_map,self.title,ax)

    #plot mean of rmse
    def plotrmsemean(self,i,ax):
        imshow_map(self.meanrmse_map,"RMSE of "+self.title,ax,vmin=0,vmax=1.5)

#世界地図とイメージの描画
def imshow_map(data,title,ax1,vmin=None,vmax=None):
    ax1.coastlines(resolution='110m')
    ax1.set_extent([-180,180,-90,90], ccrs.PlateCarree(central_longitude=-180.0))
    ax1.set_title(title)
            
    #緯度・経度線のグリッドの設定
    ax1.set_xticks([-180,-120,-60,0, 60, 120, 180], crs=ccrs.PlateCarree(central_longitude=-180.0))
    ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree(central_longitude=-180.0))
    lon_formatter = LongitudeFormatter(dateline_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter,)
    ax1.yaxis.set_major_formatter(lat_formatter)
    if vmin==None:
        im=plt.imshow(data,extent=[-180,180,-90,90],cmap=cm.jet  , origin='lower')#,vmin=220,vmax=280
    else:
        im=plt.imshow(data,extent=[-180,180,-90,90],cmap=cm.jet ,vmin=vmin,vmax=vmax , origin='lower')#,vmin=220,vmax=280

#観測点を描画 obs:raob,reg1,reg2,...,reg5
def plotobs(obs):
    x=[]
    y=[]
    s=nature['t'][0][0].shape
    print(s)
    with open("../obs/station_"+obs+".tbl","r") as f:
        f.readline()
        f.readline()
        for i in f.readlines():
            u,v=i.split()
            u,v=(float(u)-1)/s[1]*360-180,(float(v)-1)/s[0]*180-90 #stringのインデックスから緯度経度に変換
            print(u,v)
            x.append(u)
            y.append(v)

    plt.plot(x,y,"o",color="black",markersize=2)

#世界地図を描画
def worldmap(ax):
    ax.add_feature(cfeature.LAND)
    ax.coastlines(lw=0.5)
    ax.gridlines(linestyle='-', color='gray')

#cdで関数を終わると自動でcdが元に戻る便利な奴
def local_chdir(func):
    def _inner(*args, **kwargs):
        # 元のカレントディレクトリを変数に代入
        dir_original = os.getcwd()
        # 渡された関数実行
        ret = func(*args, **kwargs)
        # カレントディレクトリを元に戻す
        os.chdir(dir_original)
        return ret
    return _inner

#cdoを用いてgrads から nc4に変換
def convert_grds2nc4(control_file,outfilename):
    cmd=["cdo", "-f", "nc4", "import_binary", control_filename, outfilename]
    print(os.getcwd(),cmd)
    subprocess.run(cmd)

# ncファイルを読み込む。必要に応じて変換を行う
@local_chdir
def get_nc(path):
    meannc_filename="mean.nc"
    os.chdir(path)
    if not os.path.exists(meannc_filename): #file not exists?
        convert_grds2nc4(control_filename,meannc_filename)
    try:
        ncd= nc.Dataset(meannc_filename,"r")
    except:
        print("Error: while reading \""+path+meannc_filename+"\"", file=sys.stderr)
        exit(1)
    return ncd

#1時刻のrmse面積平均を求める
def rmse(t,level,nature,mean):
    tmean=np.array(mean.variables['t'][t][level])
    tnature=np.array(nature.variables['t'][t][level])
    tdiff=tmean-tnature
    latitude=mean.variables['lat']
    scos,srmse=0,0
    for lat,diffline in zip(latitude,tdiff):
        cos=math.cos(math.radians(lat)) #面積平均、高緯度ほど面積が小さい
        srmse+=sum(diffline**2)*cos
        scos+=cos
    rmse=math.sqrt(srmse/scos/len(tdiff[0]))#横の格子点で割る
    return rmse




sys_dirs=[ [("SYS20200815_RELAXATION_3LPFGM_M000040L"+ '{:04}'.format(i)+"IADP_RTPS0."+str(m)+"_woWSM_infMP_fac2.5_fgt1.00_rsmp002",str(m)+","+str(i)) for m in range(55,75+5,5)]for i in range(700,1100,100)]

def main():
    global nature,expr_results
    expr_results=[]
    axes = plt.axes()
    axes.set_ylim([0.4,1])
    plt.title("SPEEDY 40mem RMSE")
    plt.xlabel('RTPS')
    plt.ylabel('RMSE[K]')
    plt.grid(which='minor',color='gray',linestyle='--')
    nature=get_nc(nature_dir)
    print(nature)

#plot observation points
    # for idx,i in enumerate(["reg3","reg4","reg5"]):
    #     ax1=plt.subplot(2,2,idx+1,projection=ccrs.PlateCarree(central_longitude=-180.0))
    #     plt.title(i)
    #     worldmap(ax1)
    #     plotobs(i)
    # plt.show()
    # return

    for idx,j in enumerate(sys_dirs):
        for idx2,i in enumerate(j):
            print(i)
        
            #        ax1=plt.subplot(2,2,idx%7+1,projection=ccrs.PlateCarree(central_longitude=-180.0))
            expr=exp_result(i[0],i[1],nature)
            expr.get_mean()
            #            expr.get_kldv()
            #            expr.calc_kldvmean_map()
            expr.calc_rmse()

            #        expr.get_meant_map()
            #        expr.calc_rmsemean_map()
            #        expr.drawrmse()
            #        expr.plotrmsemean(idx,ax1)
            #        expr.plottmean(idx,ax1)
            expr_results.append(expr)
            #        cbar=plt.colorbar(orientation='horizontal')
            #        cbar.set_label("RMSE[K]")

    rmselist=[e.calc_meanrmse() for e in expr_results]
    xlist=[e.title for e in expr_results]
    plt.bar(xlist,rmselist)
    plt.show()

"""
    ax1=plt.subplot(2,2,level,projection=ccrs.PlateCarree(central_longitude=-180.0))
    imshow_map(expr_results[0].meanrmse_map-expr_results[1].meanrmse_map,"LETKF - LPFGM RMSE difference",ax1)
    cbar=plt.colorbar(orientation='horizontal')
    cbar.set_label("RMSE difference[K]")

    ax1=plt.subplot(2,2,4,projection=ccrs.PlateCarree(central_longitude=-180.0))
    imshow_map(expr_results[0].meankldv_map,"KLDV of LETKF",ax1)
    plotobs("raob")
    cbar=plt.colorbar(orientation='horizontal')
    cbar.set_label("KLDV")
#    plt.subplot(111)
#    cbar=plt.colorbar(orientation='horizontal')
#    cbar.set_label("Temp[K]")
#    plt.legend()
    plt.tight_layout()
    plt.show()
"""
    #plt.yticks(np.arange(0, 2+0.2, 0.2))
#    plt.show()
    #plt.savefig('out_graph.png', dpi=200, orientation='portrait', transparent=False, pad_inches=0.0)

if __name__ == '__main__':
    main()
