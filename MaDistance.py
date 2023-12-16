import gdal, ogr, os, osr
import pandas
import os
import shutil
import collections
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.spatial import distance
from scipy.stats import chi
import cmaps
import fileinput
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from math import radians, cos, sin, asin, sqrt
from matplotlib.colors import LinearSegmentedColormap

# 将气候栅格数据转换为csv格式
def TifToCsv(inPath_1, outPath_1):
	inPath_2 = os.path.join(inPath_1, "tmax_winter.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_Tmax_wt = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "tmax_spring.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_Tmax_sp = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "tmax_summer.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_Tmax_sm = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "tmax_autumn.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_Tmax_at = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "tmin_winter.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_Tmin_wt = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "tmin_spring.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_Tmin_sp = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "tmin_summer.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_Tmin_sm = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "tmin_autumn.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_Tmin_at = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "prec_winter.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_PPT_wt = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "prec_spring.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_PPT_sp = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "prec_summer.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_PPT_sm = array_1.flatten()

	inPath_2 = os.path.join(inPath_1, "prec_autumn.tif")
	raster_1 = gdal.Open(inPath_2)
	array_1 = raster_1.ReadAsArray()
	array_PPT_at = array_1.flatten()

	dataframe_1 = pandas.DataFrame({"Tmax_wt": array_Tmax_wt, "Tmax_sp": array_Tmax_sp, "Tmax_sm": array_Tmax_sm, "Tmax_at": array_Tmax_at, \
		"Tmin_wt": array_Tmin_wt, "Tmin_sp": array_Tmin_sp, "Tmin_sm": array_Tmin_sm, "Tmin_at": array_Tmin_at, \
		"PPT_wt": array_PPT_wt, "PPT_sp": array_PPT_sp, "PPT_sm": array_PPT_sm, "PPT_at": array_PPT_at})
	dataframe_1.to_csv(outPath_1)

# 数组转栅格
def ArrayToRaster(inPath_1, inPath_2, rows, columns, outPath_1):
	dataframe_1 = pandas.read_csv(inPath_1)
	dataframe_2 = pandas.read_csv(inPath_2)
	array_1 = dataframe_1["x"].tolist()
	array_2 = dataframe_2["x"].tolist()
	array_3 = numpy.full((rows * columns), -1)
	ID_1 =  0
	for index_1 in array_2:
		if array_1[ID_1] != float("inf"):
			array_3[index_1] = array_1[ID_1]
		ID_1 = ID_1 + 1
	array_3.resize((rows, columns))
	rasterOrigin = (73.4166666665653, 53.5833333332879)
	NoData_value = -1
	originX = rasterOrigin[0]
	originY = rasterOrigin[1]
	pixelWidth = 0.083333333
	pixelHeight = -0.083333333
	driver = gdal.GetDriverByName('GTiff')
	outRaster = driver.Create(outPath_1, columns, rows, 1, gdal.GDT_Float32)
	outRaster.SetGeoTransform([originX, pixelWidth, 0, originY, 0.0, pixelHeight])
	outband = outRaster.GetRasterBand(1)
	outband.WriteArray(array_3)
	outband.SetNoDataValue(NoData_value)				
	outRasterSRS = osr.SpatialReference()
	outRasterSRS.ImportFromEPSG(4326)
	outRaster.SetProjection(outRasterSRS.ExportToWkt())
	outband.FlushCache()

######################################################################################

# 将站点数据转换为csv格式
def TxtToCsv(inPath_1, outPath_1):
	for dirPath, dirname, filenames in os.walk(inPath_1):
		for filename in filenames:
			inPath_2 = os.path.join(inPath_1, filename)
			for line in fileinput.input(inPath_2, inplace = 1):
			    if not fileinput.isfirstline():
			        print(line.replace("\n", ""))
	for dirPath, dirname, filenames in os.walk(inPath_1):
		for filename in filenames:
			inPath_2 = os.path.join(inPath_1, filename)
			outPath_2 = os.path.join(outPath_1, filename.replace(".txt", ".csv", 1))
			txt_1 = numpy.loadtxt(inPath_2)  
			dataframe_1 = pandas.DataFrame(txt_1)  
			dataframe_1.columns = ["V01301", "V04001", "V04002", "V11046", "V11211_CHAR", "V10302", "V12012", "V10301", "V12011", "V13305", "V10004_701", \
			"V11291_701", "V12001_701", "V13004_701", "V13003_701", "V12012_701", "V12011_701", "V13353", "V14033", "V14032", "V11042", "V11296_CHAR", "V13052", "V13007"]
			dataframe_1.to_csv(outPath_2, index = False)  

# 计算各季节的温度和降水
def CalculateSeason(inPath_1, outPath_1, outPath_2):
	dataframe_1 = pandas.DataFrame(columns = ["V01301", "V04001", "V04002", "V11046", "V11211_CHAR", "V10302", "V12012", "V10301", "V12011", "V13305", "V10004_701", \
	"V11291_701", "V12001_701", "V13004_701", "V13003_701", "V12012_701", "V12011_701", "V13353", "V14033", "V14032", "V11042", "V11296_CHAR", "V13052", "V13007"])
	for dirPath, dirname, filenames in os.walk(inPath_1):
		for filename in filenames:
			inPath_2 = os.path.join(inPath_1, filename)
			dataframe_3 = pandas.read_csv(inPath_2)
			dataframe_1 = pandas.concat([dataframe_1, dataframe_3])
	dataframe_1.sort_values(by = ["V01301", "V04001", "V04002"], inplace = True)
	dataframe_1.to_csv(outPath_1, index = False)
	station_1 = dataframe_1["V01301"].unique()
	dataframe_2 = pandas.DataFrame(columns = ["Station_tag", "Tmax_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmin_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "PPT_wt", "PPT_sp", "PPT_sm", "PPT_at"])
	stationTag_1 = 0
	for station_2 in station_1:
		stationTag_1 = stationTag_1 + 1
		dataframe_3 = dataframe_1[dataframe_1["V01301"] == station_2]
		year_1 = dataframe_3["V04001"].unique()
		for year_2 in year_1:
			dataframe_4 = dataframe_3[dataframe_3["V04001"] == year_2]
			month_1 = dataframe_4["V04002"].tolist()
			tmax_1 = dataframe_4["V12011"].tolist()
			tmin_1 = dataframe_4["V12012"].tolist()
			prec_1 = dataframe_4["V13305"].tolist()
			winter_sum_tmax_1 = 0
			winter_sum_tmin_1 = 0
			winter_sum_prec_1 = 0				
			spring_sum_tmax_1 = 0
			spring_sum_tmin_1 = 0
			spring_sum_prec_1 = 0
			summer_sum_tmax_1 = 0
			summer_sum_tmin_1 = 0
			summer_sum_prec_1 = 0
			autumn_sum_tmax_1 = 0
			autumn_sum_tmin_1 = 0
			autumn_sum_prec_1 = 0
			winter_month_1 = 0			
			spring_month_1 = 0
			summer_month_1 = 0
			autumn_month_1 = 0										
			for index_1 in range(0, len(month_1)):
				if prec_1[index_1] != 999999:
					month_2 = month_1[index_1]
					if month_2 in [1, 2, 12]:
						winter_sum_tmax_1 = winter_sum_tmax_1 + tmax_1[index_1]
						winter_sum_tmin_1 = winter_sum_tmin_1 + tmin_1[index_1]
						winter_sum_prec_1 = winter_sum_prec_1 + prec_1[index_1]
						winter_month_1 = winter_month_1 + 1
					elif month_2 in [3, 4, 5]:
						spring_sum_tmax_1 = spring_sum_tmax_1 + tmax_1[index_1]
						spring_sum_tmin_1 = spring_sum_tmin_1 + tmin_1[index_1]
						spring_sum_prec_1 = spring_sum_prec_1 + prec_1[index_1]
						spring_month_1 = spring_month_1 + 1
					elif month_2 in [6, 7, 8]:
						summer_sum_tmax_1 = summer_sum_tmax_1 + tmax_1[index_1]
						summer_sum_tmin_1 = summer_sum_tmin_1 + tmin_1[index_1]
						summer_sum_prec_1 = summer_sum_prec_1 + prec_1[index_1]
						summer_month_1 = summer_month_1 + 1
					else:
						autumn_sum_tmax_1 = autumn_sum_tmax_1 + tmax_1[index_1]
						autumn_sum_tmin_1 = autumn_sum_tmin_1 + tmin_1[index_1]
						autumn_sum_prec_1 = autumn_sum_prec_1 + prec_1[index_1]
						autumn_month_1 = autumn_month_1 + 1
			if winter_month_1 != 0 and spring_month_1 != 0 and summer_month_1 != 0 and autumn_month_1 != 0:
				winter_mean_tmax = winter_sum_tmax_1 / winter_month_1
				winter_mean_tmin = winter_sum_tmin_1 / winter_month_1
				spring_mean_tmax = spring_sum_tmax_1 / spring_month_1
				spring_mean_tmin = spring_sum_tmin_1 / spring_month_1
				summer_mean_tmax = summer_sum_tmax_1 / summer_month_1
				summer_mean_tmin = summer_sum_tmin_1 / summer_month_1
				autumn_mean_tmax = autumn_sum_tmax_1 / autumn_month_1
				autumn_mean_tmin = autumn_sum_tmin_1 / autumn_month_1			
				dataframe_2 = dataframe_2.append([{"Station_tag": stationTag_1, "Tmax_wt": winter_mean_tmax, "Tmax_sp": spring_mean_tmax, "Tmax_sm": summer_mean_tmax, "Tmax_at": autumn_mean_tmax, \
					"Tmin_wt": winter_mean_tmin, "Tmin_sp": spring_mean_tmin, "Tmin_sm": summer_mean_tmin, "Tmin_at": autumn_mean_tmin, \
					"PPT_wt": winter_sum_prec_1, "PPT_sp": spring_sum_prec_1, "PPT_sm": summer_sum_prec_1, "PPT_at": autumn_sum_prec_1}])
		print(str(station_2) + " is done!")
	dataframe_2.to_csv(outPath_2, index = False)

# 创建最邻近点ID
def CreateNearestID(inPath_1, grids, outPath_1):
	dataframe_1 = pandas.read_csv(inPath_1)
	dataframe_2 = pandas.DataFrame(columns = ["V1", "V2", "V3", "V4"])
	for index_1 in range(0, grids):
		dataframe_3 = dataframe_1[dataframe_1["INPUT_FID"] == index_1]
		dataframe_3 = dataframe_3.sort_values(by = ["DISTANCE"])
		NEAR_FID = dataframe_3["NEAR_FID"].tolist()
		try:
			dataframe_2 = dataframe_2.append([{"V1": NEAR_FID[0] + 1, "V2": NEAR_FID[1] + 1, "V3": NEAR_FID[2] + 1, "V4": NEAR_FID[3] + 1}])
			print(str(index_1) + " is done!")
		except Exception as e:
			print("missing： " + str(index_1))
	dataframe_2.to_csv(outPath_1, index = False)

# 创建陆地区域ID
def CreateLandID(inPath_1, outPath_1):
	raster_1 = gdal.Open(inPath_1)
	array_1 = raster_1.ReadAsArray()
	array_2 = array_1.flatten()
	ID_1 = []
	for index_1 in range(0, len(array_2)):
		if array_2[index_1] == 1:
			ID_1.append(index_1)
	dataframe_1 = pandas.DataFrame(columns = ["x"])
	dataframe_1["x"] = ID_1
	dataframe_1.to_csv(outPath_1, index = False)

######################################################################################

# 计算马氏距离（climatic novelty）
def MaDistance1(inPath_1, inPath_2, inPath_3, inPath_4, inPath_5, outPath_1):
	dataframe_1 = pandas.read_csv(inPath_1)
	dataframe_2 = pandas.read_csv(inPath_2)
	dataframe_3 = pandas.read_csv(inPath_3)
	dataframe_4 = pandas.read_csv(inPath_4)
	dataframe_5 = pandas.read_csv(inPath_5)
	array_1 = dataframe_1.values[:, 1: 13]
	array_2 = dataframe_2.values[:, 1: 13]
	array_3 = dataframe_5.values
	ID_1 = dataframe_2.values[:, 0: 1]
	location_1 = array_3[:, 0]
	index_1 = 0
	dataframe_6 = pandas.DataFrame(columns = ["x"])
	for index_2 in location_1:
		try:
			dataframe_7 = dataframe_1[dataframe_1["x"] == index_2]
			array_4 = dataframe_7.values[:, 1: 13][0]
			dataframe_8 = dataframe_4[dataframe_4["x"] == index_2]
			ID_2 = dataframe_8.values[:, 1: 2][0][0]
			dataframe_9 = dataframe_3[dataframe_3["Station_tag"] == ID_2]
			array_5 = dataframe_9.values[:, 1: 13]
			array_6 = numpy.cov(array_5, rowvar = False)
			array_7 = numpy.linalg.inv(array_6)
			distance_1 = []
			for index_3 in range(0, len(dataframe_2)):
				distance_2 = distance.mahalanobis(array_4, array_2[index_3], array_7)
				possibility_1 = chi.cdf(distance_2, 12)
				percentile_1 = chi.ppf(possibility_1, 1)
				distance_1.append(percentile_1)
			minValue = min(distance_1)
			dataframe_6 = dataframe_6.append([{"x": minValue}])
			print(str(index_2) + " is done!")
		except Exception as e:
			print(str(index_2) + " doesn't exist!")
		index_1 = index_1 + 1
	outPath_3 = os.path.join(outPath_1, "Land_MaDistance_1.csv")
	dataframe_6.to_csv(outPath_3, index = False)

# 绘制栅格（climatic novelty）
def drawRaster1(inPath_1, rows, columns, iCmap, iTitle, outPath_1):
	data1 = gdal.Open(inPath_1)
	data2 = data1.ReadAsArray()
	data3 = numpy.full((rows, columns), numpy.nan)
	for row in range(0, rows):
		for column in range(0, columns):
			if data2[rows - 1 - row, column] != -1:
				data3[row, column] = data2[rows - 1 - row, column]
	fig = plt.figure(figsize = (8, 6.5), dpi = 330)
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 13, urcrnrlat = 49, llcrnrlon = 81, urcrnrlon = 149, resolution = "h")
	map.drawparallels(numpy.arange(-15., 85., 15.), labels = [1, 0, 0, 0], color = "grey", fontsize = 8)
	map.drawmeridians(numpy.arange(10., 180., 15.), labels = [0, 0, 0, 1], color = "grey", fontsize = 8)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	x1 = numpy.linspace(73.4166666665653, 135.166666419565, columns)
	y1 = numpy.linspace(3.8333335322879, 53.5833333332879, rows)
	x2, y2 = numpy.meshgrid(x1, y1)
	x3, y3 = map(x2, y2)
	data4 = map.pcolor(x3, y3, numpy.squeeze(data3), cmap = iCmap)
	clobar = map.colorbar(data4, location = "bottom", size = "2%", pad = "4%", ticks = [0, 0.9945, 1.96, 2.9677, 4.0128, 4.9912, 5.9978, 6.9923, 7.9868])
	clobar.set_ticklabels(["0σ", "1σ", "2σ", "3σ", "4σ", "5σ", "6σ", "7σ", "8σ"])
	clobar.ax.tick_params(labelsize = 8)
	plt.subplots_adjust(left = 0.05, bottom = 0.01, right = 0.95, top = 0.99)
	# plt.legend(loc = "lower left")
	plt.title(iTitle, fontsize = 12)
	fig.add_axes([0.7809, 0.11, 0.1761, 0.3])
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 2, urcrnrlat = 25, llcrnrlon = 107, urcrnrlon = 125, resolution = "h")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	plt.savefig(outPath_1, dpi = 330)

######################################################################################

# 计算马氏距离（climatic analogs）
def MaDistance2(inPath_1, inPath_2, inPath_3, inPath_4, inPath_5, outPath_1):
	dataframe_1 = pandas.read_csv(inPath_1)
	dataframe_2 = pandas.read_csv(inPath_2)
	dataframe_3 = pandas.read_csv(inPath_3)
	dataframe_4 = pandas.read_csv(inPath_4)
	dataframe_5 = pandas.read_csv(inPath_5)
	array_1 = dataframe_1.values[:, 1: 13]
	array_2 = dataframe_2.values[:, 1: 13]
	array_3 = dataframe_5.values
	ID_1 = dataframe_2.values[:, 0: 1]
	name_1 = array_3[:, 0]
	location_1 = array_3[:, 5]
	index_1 = 0
	dataframe_6 = pandas.DataFrame(columns = ["city", "furLon", "furLat", "conLon", "conLat", "minimum distance"])
	for index_2 in location_1:
		try:
			dataframe_7 = dataframe_1[dataframe_1["x"] == index_2]
			array_4 = dataframe_7.values[:, 1: 13][0]
			dataframe_8 = dataframe_4[dataframe_4["x"] == index_2]
			ID_2 = dataframe_8.values[:, 1: 2][0][0]
			dataframe_9 = dataframe_3[dataframe_3["Station_tag"] == ID_2]
			array_5 = dataframe_9.values[:, 1: 13]
			array_6 = numpy.cov(array_5, rowvar = False)
			array_7 = numpy.linalg.inv(array_6)
			dataframe_10 = pandas.DataFrame(columns = ["x"])
			distance_1 = []
			for index_3 in range(0, len(dataframe_2)):
				distance_2 = distance.mahalanobis(array_4, array_2[index_3], array_7)
				possibility_1 = chi.cdf(distance_2, 12)
				percentile_1 = chi.ppf(possibility_1, 1)
				distance_1.append(percentile_1)
			minValue = min(distance_1)
			minIndex = distance_1.index(min(distance_1))
			minLocation = ID_1[minIndex][0]
			furRow = int(index_2 / 741) + 1
			furColumn = int(index_2 % 741)
			furLon = furColumn * 0.083333333 + 73.4166666665653
			furLat = 53.5833333332879 - furRow * 0.083333333
			conRow = int(minLocation / 741) + 1
			conColumn = int(minLocation % 741)
			conLon = conColumn * 0.083333333 + 73.4166666665653
			conLat = 53.5833333332879 - conRow * 0.083333333
			dataframe_6 = dataframe_6.append([{"city": name_1[index_1], "furLon": furLon, "furLat": furLat, "conLon": conLon, "conLat": conLat, "minimum distance": minValue}])
			dataframe_10["x"] = distance_1
			outPath_2 = os.path.join(outPath_1, "MaDistance_" + name_1[index_1] + ".csv")
			dataframe_10.to_csv(outPath_2, index = False)
			print(str(index_2) + " is done!")
		except Exception as e:
			print(str(index_2) + " doesn't exist!")
		index_1 = index_1 + 1
	outPath_3 = os.path.join(outPath_1, "All_MaDistance_1.csv")
	dataframe_6.to_csv(outPath_3, index = False)

# 绘制栅格（climatic analogs）
def drawRaster2(inPath_1, inPath_2, city, rows, columns, iCmap, iTitle, outPath_1):
	data1 = gdal.Open(inPath_1)
	data2 = data1.ReadAsArray()
	data3 = numpy.full((rows, columns), numpy.nan)
	for row in range(0, rows):
		for column in range(0, columns):
			if data2[rows - 1 - row, column] != -1:
				data3[row, column] = data2[rows - 1 - row, column]
	fig = plt.figure(figsize = (8, 6.5), dpi = 330)
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 13, urcrnrlat = 49, llcrnrlon = 81, urcrnrlon = 149, resolution = "h")
	map.drawparallels(numpy.arange(-15., 85., 15.), labels = [1, 0, 0, 0], color = "grey", fontsize = 8)
	map.drawmeridians(numpy.arange(10., 180., 15.), labels = [0, 0, 0, 1], color = "grey", fontsize = 8)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	dataframe_1 = pandas.read_csv(inPath_2)
	dataframe_2 = dataframe_1[dataframe_1["city"] == city]
	array_1 = dataframe_2.values[0]
	map.plot(array_1[1], array_1[2], ".", markersize = 2, color = "purple", latlon = True)
	# map.plot(array_1[3], array_1[4], "bo", markersize = 1, color = "b", latlon = True, label = "climatic analog")
	x_1, y_1 = map(array_1[1], array_1[2])
	x_2, y_2 = map(array_1[3], array_1[4])
	plt.arrow(x_1, y_1, x_2 - x_1, y_2 - y_1, linewidth = 0.6, length_includes_head = True, color = "purple", head_width = 30000, head_length = 50000)
	map.plot(array_1[3], array_1[4], ".", markersize = 2, color = "purple", latlon = True)
	# x_3, y_3 = map(array_1[1] - 5, array_1[2] + 1)
	x_3, y_3 = map(array_1[1] + 1, array_1[2] + 1)
	# x_3, y_3 = map(array_1[1] - 5, array_1[2] + 1)
	# x_3, y_3 = map(array_1[1] + 1, array_1[2])
	plt.text(x_3, y_3, city + " (" + str(round(array_1[5], 2)) + ")", fontsize = 10, color = "purple")
	# x_4, y_4 = map(array_1[3] - 2, array_1[4] - 1)
	x_4, y_4 = map(array_1[3] + 1, array_1[4] - 1)
	# x_4, y_4 = map(array_1[3] - 3, array_1[4] - 1)
	# x_4, y_4 = map(array_1[3] + 1, array_1[4])
	plt.text(x_4, y_4, array_1[6], fontsize = 10, color = "purple")
	x1 = numpy.linspace(73.4166666665653, 135.166666419565, columns)
	y1 = numpy.linspace(3.8333335322879, 53.5833333332879, rows)
	x2, y2 = numpy.meshgrid(x1, y1)
	x3, y3 = map(x2, y2)
	data4 = map.pcolor(x3, y3, numpy.squeeze(data3), cmap = iCmap)
	clobar = map.colorbar(data4, location = "bottom", size = "2%", pad = "4%", ticks = [0, 0.9945, 1.96, 2.9677, 4.0128, 4.9912, 5.9978, 6.9923, 7.9868])
	clobar.set_ticklabels(["0σ", "1σ", "2σ", "3σ", "4σ", "5σ", "6σ", "7σ", "8σ"])
	clobar.ax.tick_params(labelsize = 8)
	plt.subplots_adjust(left = 0.05, bottom = 0.01, right = 0.95, top = 0.99)
	# plt.legend(loc = "lower left")
	plt.title(iTitle, fontsize = 12)
	fig.add_axes([0.7809, 0.11, 0.1761, 0.3])
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 2, urcrnrlat = 25, llcrnrlon = 107, urcrnrlon = 125, resolution = "h")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	plt.savefig(outPath_1, dpi = 330)

# 绘制类似气候（climatic analogs）
def drawLink(inPath_1, iTitle, outPath_1):
	fig = plt.figure(figsize = (8, 6.5), dpi = 330)
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 13, urcrnrlat = 49, llcrnrlon = 81, urcrnrlon = 149, resolution = "h")
	map.drawparallels(numpy.arange(-15., 85., 15.), labels = [1, 0, 0, 0], color = "grey", fontsize = 8)
	map.drawmeridians(numpy.arange(10., 180., 15.), labels = [0, 0, 0, 1], color = "grey", fontsize = 8)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	dataframe_1 = pandas.read_csv(inPath_1)
	array_1 = dataframe_1.values
	# map.plot(array_1[:, 1], array_1[:, 2], "o", markersize = 0.5, color = "k", latlon = True)
	# map.plot(array_1[:, 3], array_1[:, 4], "bo", markersize = 1, color = "y", latlon = True, label = "climatic analog")
	eastTag = 0
	southTag = 0
	westTag = 0
	northTag = 0
	distance_1 = 0
	for index_1 in range(0, len(array_1)):
		x_1, y_1 = map(array_1[index_1][1], array_1[index_1][2])
		x_2, y_2 = map(array_1[index_1][3], array_1[index_1][4])
		lon1 = radians(array_1[index_1][1])
		lat1 = radians(array_1[index_1][2])
		lon2 = radians(array_1[index_1][3])
		lat2 = radians(array_1[index_1][4])
		dlon = lon2 - lon1 
		dlat = lat2 - lat1 
		a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
		c = 2 * asin(sqrt(a)) 
		r = 6371 # 地球平均半径，单位为公里
		distance_1 = distance_1 + c * r
		if array_1[index_1][1] > array_1[index_1][3]:
			if array_1[index_1][2] > array_1[index_1][4]:
				map.plot(array_1[index_1, 1], array_1[index_1, 2], "o", markersize = 0.5, color = "b", latlon = True)
				p_1 = plt.arrow(x_1, y_1, x_2 - x_1, y_2 - y_1, linewidth = 0.5, length_includes_head = True, color = "b", head_width = 30000, head_length = 30000)
				southTag = southTag + 1
				westTag = westTag + 1
			else:
				map.plot(array_1[index_1, 1], array_1[index_1, 2], "o", markersize = 0.5, color = "r", latlon = r)
				p_2 = plt.arrow(x_1, y_1, x_2 - x_1, y_2 - y_1, linewidth = 0.5, length_includes_head = True, color = "r", head_width = 30000, head_length = 30000)
				northTag = northTag + 1
				westTag = westTag + 1
		else:
			if array_1[index_1][2] > array_1[index_1][4]:
				map.plot(array_1[index_1, 1], array_1[index_1, 2], "o", markersize = 0.5, color = "g", latlon = True)
				p_3 = plt.arrow(x_1, y_1, x_2 - x_1, y_2 - y_1, linewidth = 0.5, length_includes_head = True, color = "g", head_width = 30000, head_length = 30000)
				southTag = southTag + 1
				eastTag = eastTag + 1
			else:
				map.plot(array_1[index_1, 1], array_1[index_1, 2], "o", markersize = 0.5, color = "y", latlon = True)
				p_4 = plt.arrow(x_1, y_1, x_2 - x_1, y_2 - y_1, linewidth = 0.5, length_includes_head = True, color = "y", head_width = 30000, head_length = 30000)	
				northTag = northTag + 1
				eastTag = eastTag + 1
	legend_elements = [Line2D([0], [0], color = "b", linewidth = 6, label = "southwest"), Line2D([0], [0], color = "r", linewidth = 6, label = "northwest"), \
	Line2D([0], [0], color = "g", linewidth = 6, label = "southeast"), Line2D([0], [0], color = "y", linewidth = 6, label = "northeast")]
	plt.legend(handles = legend_elements, loc = "lower left")
	plt.subplots_adjust(left = 0.05, bottom = 0.01, right = 0.95, top = 0.99)
	plt.title(iTitle, fontsize = 12)
	fig.add_axes([0.7809, 0.085, 0.1761, 0.3])
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 2, urcrnrlat = 25, llcrnrlon = 107, urcrnrlon = 125, resolution = "h")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	plt.savefig(outPath_1, dpi = 330)

# 绘制类似气候（climatic analogs）
def drawDis(inPath_1, iTitle, outPath_1):
	fig = plt.figure(figsize = (8, 6.5), dpi = 330)
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 13, urcrnrlat = 49, llcrnrlon = 81, urcrnrlon = 149, resolution = "h")
	map.drawparallels(numpy.arange(-15., 85., 15.), labels = [1, 0, 0, 0], color = "grey", fontsize = 8)
	map.drawmeridians(numpy.arange(10., 180., 15.), labels = [0, 0, 0, 1], color = "grey", fontsize = 8)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	dataframe_1 = pandas.read_csv(inPath_1)
	array_1 = dataframe_1.values
	map.plot(array_1[:, 3], array_1[:, 4], "bo", markersize = 1, color = "y", latlon = True, label = "climatic analog")
	for index_1 in range(0, len(array_1)):
		x_1, y_1 = map(array_1[index_1][1], array_1[index_1][2])
		x_2, y_2 = map(array_1[index_1][3], array_1[index_1][4])
		if array_1[index_1][5] <= 1.96:
			if array_1[index_1][1] > array_1[index_1][3]:
				if array_1[index_1][2] > array_1[index_1][4]:
					map.plot(array_1[index_1][1], array_1[index_1][2], "v", markersize = array_1[index_1][5], color = "b", latlon = True)
				else:
					map.plot(array_1[index_1][1], array_1[index_1][2], "v", markersize = array_1[index_1][5], color = "r", latlon = True)
			else:
				if array_1[index_1][2] > array_1[index_1][4]:
					map.plot(array_1[index_1][1], array_1[index_1][2], "v", markersize = array_1[index_1][5], color = "g", latlon = True)
				else:
					map.plot(array_1[index_1][1], array_1[index_1][2], "v", markersize = array_1[index_1][5], color = "y", latlon = True)
		else:
			if array_1[index_1][1] > array_1[index_1][3]:
				if array_1[index_1][2] > array_1[index_1][4]:
					map.plot(array_1[index_1][1], array_1[index_1][2], "o", markersize = array_1[index_1][5], color = "b", latlon = True)
				else:
					map.plot(array_1[index_1][1], array_1[index_1][2], "o", markersize = array_1[index_1][5], color = "r", latlon = True)
			else:
				if array_1[index_1][2] > array_1[index_1][4]:
					map.plot(array_1[index_1][1], array_1[index_1][2], "o", markersize = array_1[index_1][5], color = "g", latlon = True)
				else:
					map.plot(array_1[index_1][1], array_1[index_1][2], "o", markersize = array_1[index_1][5], color = "y", latlon = True)
	legend_elements = [Line2D([0], [0], color = "b", linewidth = 6, label = "southwest"), Line2D([0], [0], color = "r", linewidth = 6, label = "northwest"), \
	Line2D([0], [0], color = "g", linewidth = 6, label = "southeast"), Line2D([0], [0], color = "y", linewidth = 6, label = "northeast"), \
	Line2D([0], [0], marker = "v", color = "w", label = "≤ 2σ", markerfacecolor = "k", markersize = 6), Line2D([0], [0], marker = "o", color = "w", label = "> 2σ", markerfacecolor = "k", markersize = 6)]
	plt.legend(handles = legend_elements, loc = "lower left")
	plt.subplots_adjust(left = 0.05, bottom = 0.01, right = 0.95, top = 0.99)
	plt.title(iTitle, fontsize = 12)
	fig.add_axes([0.7809, 0.085, 0.1761, 0.3])
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 2, urcrnrlat = 25, llcrnrlon = 107, urcrnrlon = 125, resolution = "h")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	plt.savefig(outPath_1, dpi = 330)
######################################################################################

# 计算马氏距离（livable area）
def MaDistance3(inPath_1, inPath_2, inPath_3, inPath_4, inPath_5, outPath_1):
	dataframe_1 = pandas.read_csv(inPath_1) # 当代气候
	dataframe_2 = pandas.read_csv(inPath_2) # 未来气候
	dataframe_3 = pandas.read_csv(inPath_3) # 气象站点
	dataframe_4 = pandas.read_csv(inPath_4) # 最近气象站点
	dataframe_5 = pandas.read_csv(inPath_5) # 计算区域
	array_2 = dataframe_2.values[:, 1: 13]
	array_3 = dataframe_5.values
	location_1 = array_3[:, 0]
	dataframe_6 = pandas.DataFrame(columns = ["x"])
	distance_1 = numpy.full((len(dataframe_2), len(location_1)), float("inf"))
	for index_1 in range(0, len(location_1)):
		try:
			dataframe_7 = dataframe_1[dataframe_1["x"] == location_1[index_1]]
			array_4 = dataframe_7.values[:, 1: 13][0]
			dataframe_8 = dataframe_4[dataframe_4["x"] == location_1[index_1]]
			ID_2 = dataframe_8.values[:, 1: 2][0][0]
			dataframe_9 = dataframe_3[dataframe_3["Station_tag"] == ID_2]
			array_5 = dataframe_9.values[:, 1: 13]
			array_6 = numpy.cov(array_5, rowvar = False)
			array_7 = numpy.linalg.inv(array_6)
			for index_2 in range(0, len(dataframe_2)):
				distance_2 = distance.mahalanobis(array_4, array_2[index_2], array_7)
				distance_1[index_2][index_1] = distance_2
			print(str(index_1) + " is done!")
		except Exception as e:
			print(str(index_1) + " doesn't exist!")
	distance_3 = []
	for index_4 in range(0, len(dataframe_2)):
		distance_4 = distance_1[index_4].tolist()
		possibility_1 = chi.cdf(min(distance_4), 12)
		percentile_1 = chi.ppf(possibility_1, 1)
		distance_3.append(percentile_1) 
	dataframe_6["x"] = distance_3
	dataframe_6.to_csv(outPath_1, index = False)

# 绘制栅格（livable area）
def drawRaster3(inPath_1, rows, columns, iCmap, iTitle, outPath_1):
	data1 = gdal.Open(inPath_1)
	data2 = data1.ReadAsArray()
	data3 = numpy.full((rows, columns), numpy.nan)
	for row in range(0, rows):
		for column in range(0, columns):
			if data2[rows - 1 - row, column] != 0:
				data3[row, column] = 1
	fig = plt.figure(figsize = (8, 6.5), dpi = 330)
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 13, urcrnrlat = 49, llcrnrlon = 81, urcrnrlon = 149, resolution = "h")
	map.drawparallels(numpy.arange(-15., 85., 15.), labels = [1, 0, 0, 0], color = "grey", fontsize = 8)
	map.drawmeridians(numpy.arange(10., 180., 15.), labels = [0, 0, 0, 1], color = "grey", fontsize = 8)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	x1 = numpy.linspace(73.4166666665653, 135.166666419565, columns)
	y1 = numpy.linspace(3.8333335322879, 53.5833333332879, rows)
	x2, y2 = numpy.meshgrid(x1, y1)
	x3, y3 = map(x2, y2)
	data4 = map.pcolor(x3, y3, numpy.squeeze(data3), cmap = iCmap, label = "residential area")
	legend_elements = [Line2D([0], [0], color = "yellow", linewidth = 6, label = "residential area")]
	plt.legend(handles = legend_elements, loc = "lower left")
	plt.subplots_adjust(left = 0.05, bottom = 0.01, right = 0.95, top = 0.99)
	plt.title(iTitle, fontsize = 12)
	fig.add_axes([0.7809, 0.085, 0.1761, 0.3])
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 2, urcrnrlat = 25, llcrnrlon = 107, urcrnrlon = 125, resolution = "h")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	plt.savefig(outPath_1, dpi = 330)

# 绘制栅格（climatic novelty）
def drawRaster4(inPath_1, rows, columns, iCmap, iTitle, outPath_1):
	data1 = gdal.Open(inPath_1)
	data2 = data1.ReadAsArray()
	data3 = numpy.full((rows, columns), numpy.nan)
	for row in range(0, rows):
		for column in range(0, columns):
			if data2[rows - 1 - row, column] != -1:
				data3[row, column] = data2[rows - 1 - row, column]
	fig = plt.figure(figsize = (8, 6.5), dpi = 330)
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 13, urcrnrlat = 49, llcrnrlon = 81, urcrnrlon = 149, resolution = "h")
	map.drawparallels(numpy.arange(-15., 85., 15.), labels = [1, 0, 0, 0], color = "grey", fontsize = 8)
	map.drawmeridians(numpy.arange(10., 180., 15.), labels = [0, 0, 0, 1], color = "grey", fontsize = 8)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	x1 = numpy.linspace(73.4166666665653, 135.166666419565, columns)
	y1 = numpy.linspace(3.8333335322879, 53.5833333332879, rows)
	x2, y2 = numpy.meshgrid(x1, y1)
	x3, y3 = map(x2, y2)
	data4 = map.pcolor(x3, y3, numpy.squeeze(data3), cmap = iCmap)
	clobar = map.colorbar(data4, location = "bottom", size = "2%", pad = "4%", ticks = [0, 1.96, 4.0128, 7.9868])
	clobar.set_ticklabels(["Highly similar", "Moderately similar", "Lightly similar", "Dissimilar"])
	clobar.ax.tick_params(labelsize = 8)
	plt.subplots_adjust(left = 0.05, bottom = 0.01, right = 0.95, top = 0.99)
	plt.title(iTitle, fontsize = 12)
	fig.add_axes([0.7809, 0.11, 0.1761, 0.3])
	map = Basemap(projection = "aea", lat_1 = 25, lat_2 = 47, lon_0 = 105, llcrnrlat = 2, urcrnrlat = 25, llcrnrlon = 107, urcrnrlon = 125, resolution = "h")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_province", "China_province", linewidth = 0.3, color = "grey")
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_country", "China_country", linewidth = 0.5)
	map.readshapefile(r"I:\3_climatic_analog\1_data\1_source\0_base\1_boundary\China_nine_dotted_line", "China_nine_dotted_line", linewidth = 0.5, color = "grey")	
	map.drawlsmask()
	plt.savefig(outPath_1, dpi = 330)
######################################################################################

def main():

	# # 将站点数据转换为csv格式
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\3_observation record\2_delete_row"
	# outPath_1 = r"I:\3_climatic_analog\1_data\1_source\3_observation record\3_csv"
	# TxtToCsv(inPath_1, outPath_1)

	# # 计算各季节的温度和降水
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\3_observation record\3_csv"
	# outPath_1 = r"I:\3_climatic_analog\1_data\1_source\3_observation record\4_merge\observation record_merge.csv"
	# outPath_2 = r"I:\3_climatic_analog\1_data\1_source\4_input\1_observation record.csv"
	# CalculateSeason(inPath_1, outPath_1, outPath_2)

	# # 创建最邻近点ID（climatic novelty, climatic analogs）
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\0_base\4_land_f-c\point_distance_f-c.csv"
	# outPath_1 = r"I:\3_climatic_analog\1_data\1_source\4_input\2_grid nearest ID_f-c.csv"
	# grids = 138341
	# CreateNearestID(inPath_1, grids, outPath_1)

	# # 创建最邻近点ID（livable area）
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\0_base\5_land_c-f\point_distance_c-f.csv"
	# outPath_1 = r"I:\3_climatic_analog\1_data\1_source\4_input\3_grid nearest ID_c-f.csv"
	# grids = 656
	# CreateNearestID(inPath_1, grids, outPath_1)

	# # 创建陆地区域ID（climatic novelty, climatic analogs）
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\0_base\4_land_f-c\land.tif"
	# outPath_1 = r"I:\3_climatic_analog\1_data\1_source\4_input\4_land ID_f-c.csv"
	# CreateLandID(inPath_1, outPath_1)

	# # 创建陆地区域ID（livable area）
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\0_base\5_land_c-f\China_population_density_4.tif"
	# outPath_1 = r"I:\3_climatic_analog\1_data\1_source\4_input\5_land ID_c-f.csv"
	# CreateLandID(inPath_1, outPath_1)

######################################################################################	

	# climatic novelty

	# # 将当代气候栅格数据转为csv
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\1_contemporary climate\1_1970-2000\2_5'\1_season_tif"
	# outPath_1 = r"I:\3_climatic_analog\1_data\1_source\1_contemporary climate\1_1970-2000\2_5'\2_year_csv\contemporary climate.csv"
	# TifToCsv(inPath_1, outPath_1)

	# RCP = "RCP 4.5"
	# RCPPath = "1_RCP 4.5"

	# # 将未来气候栅格数据转为csv
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\1_source\2_future climate", RCPPath, r"4_season_tif\2050")
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\1_source\2_future climate", RCPPath, r"5_year_csv\future climate_2050.csv")
	# TifToCsv(inPath_1, outPath_1)

	# # 计算马氏距离
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\1_source\2_future climate", RCPPath, r"5_year_csv\future climate_land_2080.csv")
	# inPath_2 = r"I:\3_climatic_analog\1_data\1_source\1_contemporary climate\1_1970-2000\3_1°\2_year_csv\contemporary climate_land.csv"
	# inPath_3 = r"I:\3_climatic_analog\1_data\1_source\4_input\1_observation record.csv"
	# inPath_4 = r"I:\3_climatic_analog\1_data\1_source\4_input\2_grid nearest ID_f-c.csv"
	# inPath_5 = r"I:\3_climatic_analog\1_data\1_source\4_input\4_land ID_f-c.csv"
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\1_climatic novelty", RCPPath, "1_csv")
	# MaDistance1(inPath_1, inPath_2, inPath_3, inPath_4, inPath_5, outPath_1)

	# # 数组转栅格
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\1_climatic novelty", RCPPath, r"1_csv\Land_MaDistance_1.csv")
	# inPath_2 = r"I:\3_climatic_analog\1_data\1_source\4_input\4_land ID_f-c.csv"
	# rows = 597
	# columns = 741
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\1_climatic novelty", RCPPath, r"2_tif\Land_MaDistance_1.tif")
	# ArrayToRaster(inPath_1, inPath_2, rows, columns, outPath_1)

	# 绘制栅格
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\1_climatic novelty", RCPPath, r"2_tif\Land_MaDistance_1.tif")
	# rows = 597
	# columns = 741
	# iCmap = "YlOrRd"
	# iTitle = "Distribution of climatic novelty (2080s, " + RCP + ")"
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\5_figure\2", "1_" + iTitle + ".png")
	# drawRaster1(inPath_1, rows, columns, iCmap, iTitle, outPath_1)

######################################################################################		

	# climatic analogs

	# RCP = "RCP 4.5"
	# RCPPath = "1_RCP 4.5"
	# city = "Shanghai"

	# # 计算马氏距离
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\1_source\2_future climate", RCPPath, r"5_year_csv\future climate_land_2080.csv")
	# inPath_2 = r"I:\3_climatic_analog\1_data\1_source\1_contemporary climate\1_1970-2000\3_1°\2_year_csv\contemporary climate_land.csv"
	# inPath_3 = r"I:\3_climatic_analog\1_data\1_source\4_input\1_observation record.csv"
	# inPath_4 = r"I:\3_climatic_analog\1_data\1_source\4_input\2_grid nearest ID_f-c.csv"
	# inPath_5 = r"I:\3_climatic_analog\1_data\1_source\0_base\2_cities\China_cities.csv"
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\2_climatic analogs", RCPPath, "1_csv")
	# MaDistance2(inPath_1, inPath_2, inPath_3, inPath_4, inPath_5, outPath_1)

	# # 数组转栅格
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\2_climatic analogs", RCPPath, "1_csv", "MaDistance_" + city + ".csv")
	# inPath_2 = r"I:\3_climatic_analog\1_data\1_source\4_input\4_land ID_f-c.csv"
	# rows = 597
	# columns = 741
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\2_climatic analogs", RCPPath, "2_tif", "MaDistance_" + city + ".tif")
	# ArrayToRaster(inPath_1, inPath_2, rows, columns, outPath_1)

	# 绘制栅格
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\2_climatic analogs", RCPPath, "2_tif", "MaDistance_" + city + ".tif")
	# inPath_2 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\2_climatic analogs", RCPPath, "1_csv", "All_MaDistance_3.csv")
	# rows = 597
	# columns = 741
	# iCmap = "PuBuGn"
	# iTitle = "Climatic analog maps for " + city + " (2080s, " + RCP + ")"
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\5_figure\2", "4_" + iTitle + ".png")
	# drawRaster2(inPath_1, inPath_2, city, rows, columns, iCmap, iTitle, outPath_1)

	# # 绘制气候类似
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\2_climatic analogs", RCPPath, "1_csv", "All_MaDistance_3.csv")
	# iTitle = "Distance and direction to the best climatic analog (2080s, " + RCP + ")"
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\5_figure\2", "8_" + iTitle + ".png")
	# drawLink(inPath_1, iTitle, outPath_1)

	# 绘制不相似度
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\2_climatic analogs", RCPPath, "1_csv", "All_MaDistance_3.csv")
	# iTitle = "Strength of analogy for contemporary climatic analogs (2080s, " + RCP + ")"
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\5_figure\2", "9_" + iTitle + ".png")
	# drawDis(inPath_1, iTitle, outPath_1)

######################################################################################		

	# livable area

	RCP = "RCP 4.5"
	RCPPath = "1_RCP 4.5"

	# # 计算马氏距离
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\1_contemporary climate\1_1970-2000\3_1°\2_year_csv\contemporary climate_land.csv"
	# inPath_2 = os.path.join(r"I:\3_climatic_analog\1_data\1_source\2_future climate", RCPPath, r"5_year_csv\future climate_land_2050.csv")
	# inPath_3 = r"I:\3_climatic_analog\1_data\1_source\4_input\1_observation record.csv"
	# inPath_4 = r"I:\3_climatic_analog\1_data\1_source\4_input\3_grid nearest ID_c-f.csv"
	# inPath_5 = r"I:\3_climatic_analog\1_data\1_source\4_input\5_land ID_c-f.csv"
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\3_livable area", RCPPath, "1_csv", "Livable_MaDistance_2050.csv")
	# MaDistance3(inPath_1, inPath_2, inPath_3, inPath_4, inPath_5, outPath_1)

	# # 数组转栅格
	# inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\3_livable area", RCPPath, "1_csv", "Livable_MaDistance_2050.csv")
	# inPath_2 = r"I:\3_climatic_analog\1_data\1_source\4_input\4_land ID_f-c.csv"
	# rows = 597
	# columns = 741
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\3_livable area", RCPPath, "2_tif", "Livable_MaDistance_2050.tif")
	# ArrayToRaster(inPath_1, inPath_2, rows, columns, outPath_1)

	# # 绘制栅格
	# inPath_1 = r"I:\3_climatic_analog\1_data\1_source\0_base\5_land_c-f\China_population_density_2.tif"
	# rows = 597
	# columns = 741
	# iCmap = "inferno_r"
	# iTitle = "Contemporary residential area"
	# outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\5_figure\2", "11_" + iTitle + ".png")
	# drawRaster3(inPath_1, rows, columns, iCmap, iTitle, outPath_1)

	# 绘制栅格
	inPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\3_livable area", RCPPath, "2_tif", "Livable_MaDistance_2080.tif")
	rows = 597
	columns = 741         
	iCmap = "inferno_r"
	iTitle = "Prediction of similar climate area (2080s, " + RCP + ")"
	outPath_1 = os.path.join(r"I:\3_climatic_analog\1_data\2_result\5_figure\2", "17_" + iTitle + ".png")
	drawRaster4(inPath_1, rows, columns, iCmap, iTitle, outPath_1)

main()
