#ifndef DECODE_H
#define DECODE_H
#include "Global.h"

//获取一帧中坐标为(i,j)的宏块数据
int GetDecodeBlockData(int x,int y,int16 encodeLumaData[XX][YY]);

//对I帧数据进行解码。	
int DecodeIFrame(int16 encodeLumaData[XX][YY],FrameParament &curFrameData);

//对P帧数据进行编码。
int DecodePFrame(int index,FrameParament &curFrameData,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS]);

//先进行LDPC解码然后通过Hash校验得出是否是合适的边信息。
int SyndromeDecode(int x,int y);

//P帧中设置QuantBlockData,encodeLumaData数组;
int SetQuantBlockData(int x,int y,int16 uPrismLumaData[X][Y][UPRISMBITS],
					  int16 encodeLumaData[XX][YY]);

//对一个宏块块进行解码。
int DecodeMacroBlock();	

int IQuantAndIDCT(int x,int y,int16 source[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH]);

//解码后的宏块保存到curDecodedFrameData.LumaData[XX][YY]中;
int SetDecodeData(int x,int y,FrameParament &curFrameData);

//将一帧已经解码的数据存到DECODE_FILE中
int WriteOneDecodeFrame(int frameIndex,FILE *fp1,FrameParament curFrameData);

//获取一个伴随子宏块
int GetSyndromeBits(int x,int y,uint8 syndromeLumaBits[X][Y][rows]);

//得到最佳伴随子信息
int GetBestSideInfo(int x,int y,int16 syndromeCheck[X][Y*CHECKBITS]);

//LDPC解码
int LDPCDecode();



//先通过Hash校验得出边信息再进行LDPC解码。
int SyndromeDecode3();
int SAD(int x,int y,int dx,int dy,int16 syndromeCheck[X][Y*CHECKBITS]);
int SHD(int x,int y,int dx,int dy,int16 syndromeCheck[X][Y*CHECKBITS]);
void SaveDecodeFile(int i,int j);
//LDPC解码要用的函数
void qin2tin(double tin[rows][maxCheckDegree],double qin[cols][maxVariableDegree],double q0[cols]);
void tin2qin(double qin[cols][maxVariableDegree],double tin[rows][maxCheckDegree],uint8 *syndrome);
double Operator(double *temp, int *temp1,int j,bool flag);
void GetSyndrome(uint8 *x,uint8 *syndrome);
bool CompareSyndrome(uint8 *syndrome,uint8 *xxSyndrome);

#endif