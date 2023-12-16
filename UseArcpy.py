# -- coding: UTF-8 --

import arcpy
from arcpy import env
from arcpy.sa import *
import os

# 提取中国当代气候数据
def ExtractConChina(inPath_1, outPath_1):
	for dirPath, dirname, filenames in os.walk(inPath_1):
		for filename in filenames:
			inPath_2 = os.path.join(dirPath, filename)
			inPath_3 = r"F:\3_THU\climatic analogs\data\1_source\0_base\1_boundary\China_country.shp"
			outPath_2 = os.path.join(outPath_1, filename[-16: ])
			extract_1 = ExtractByMask(inPath_2, inPath_3)
			arcpy.Resample_management(extract_1, outPath_2, "8.33333333333333E-02", "NEAREST")
			print(outPath_2)

# 提取中国未来气候数据
def ExtractFurChina(inPath_1, outPath_1):
	walk = arcpy.da.Walk(inPath_1, type = "GRID")
	for dirPath, dirname, filenames in walk:
		for filename in filenames:
			inPath_2 = os.path.join(dirPath, filename)
			inPath_3 = r"F:\3_THU\climatic analogs\data\1_source\0_base\1_boundary\China_country.shp"
			if filename[0: 4] in ["tmax", "tmin"]:
				inPath_4 = Raster(inPath_2) / 10
			else:
				inPath_4 = inPath_2
			if len(filename) == 6:
				outPath_2 = os.path.join(outPath_1, filename[0: 4] + dirPath[-40: -35] + "-" + filename[5: 6].zfill(2) + ".tif")
			else:
				outPath_2 = os.path.join(outPath_1, filename[0: 4] + dirPath[-40: -35] + "-" + filename[5: 7] + ".tif")
			extract_1 = ExtractByMask(inPath_4, inPath_3)
			extract_1.save(outPath_2)
			print(outPath_2)

def main():
# # 提取中国当代气候数据
# 	inPath_1 = r"F:\3_THU\climatic analogs\2_data\1_contemporary climate\1_download"
# 	outPath_1 = r"F:\3_THU\climatic analogs\2_data\1_contemporary climate\2_extract_resample"
# 	ExtractConChina(inPath_1, outPath_1)

# 提取中国未来气候数据
	inPath_1 = r"C:\Users\eraer\Downloads"
	outPath_1 = r"F:\3_THU\climatic analogs\data\1_source\2_future climate\2_RCP 8.5\2_extract"
	ExtractFurChina(inPath_1, outPath_1)


main()
