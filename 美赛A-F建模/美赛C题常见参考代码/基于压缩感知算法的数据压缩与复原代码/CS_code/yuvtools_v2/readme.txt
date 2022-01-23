这个是从网上下的，一些处理视频序列的MFC小程序，挺好用的。

yuvtools_v2 What's new:
2007.09.01
＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
新增：
BMP2AVI.exe
将一副副的BMP图片合并成标准windows avi文件

BMPSeg.exe
从BMP图片中分割出一个区域，支持批量操作。

SEQ2AVI.exe
将yuv420序列转换成avi文件

SeqSad.exe
求两个序列之间,每个宏块之间的差异sad
如果没有差异就不输出。
note:该程序用于排错。

msu_vmt.exe
视频序列psnr对比分析工具


更新：
BMP2GBMP.exe
BMP2SEQ.exe
YUVviewerPlus.exe

＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝

yuvtools_v1
缩写解释：
BMP:24位的BMP图片
GBMP:一组rgb24图片（不含文件头）
SEQ:4:2:0的YUV序列
YUV:单个的4:2:0的YUV文件
DYUV:分离的4:2:0的Y，U和V文件
-------------------

YUVviewerPlus.exe
对YUVviewer进行了修改，增加一下内容：
1、增加支持的格式：yuv4:4:4, yuv4:2:2, gbmp
2、增加zoom的范围

BMP2GBMP.exe
将一副副的BMP图片合并成一个没有BMP头信息的文件
note:仅支持24位bmp图片

DYUV2SEQ.exe
实现分离的yuv文件转换成YUV序列 4：2：0

ShowDIB.exe
BMP图片显示程序，多文档框架

YUV2BMP.exe
实现了YUV转换成24位的BMP图片，实现了批量转换

BMP2SEQ.exe
将一系列24位或8位的BMP图片转换成4：2：0的YUV序列。

DYUV2BMP.exe
将分离的Y，U，V转换成24位的BMP图片，实现了批量转换

GBMP2SEQ.exe
实现包含一组rgb24图片（不含文件头）的单一文件到yuv4:2:0序列文件的转换

SEQ2BMP.exe
实现了SEQ2BMP的程序
输出BMP文件为24位真彩

SeqCut.exe
实现对YUV4:2:0文件的剪切操作
即从序列文件中取出一段序列

SeqSnr.exe
实现了两序列对应帧之间Y分量的SNR求取，并给出平均值

YUV2SEQ.exe
将单帧的YUV文件转换位YUV序列 4：2：0

YUV2SEQ2.exe
将单帧的YUV文件转换位YUV序列 4：2：0
可以选择目标图像的位置和大小


========================
Peter Lee
lspbeyond@sohu.com
updated on 2005.04.14

