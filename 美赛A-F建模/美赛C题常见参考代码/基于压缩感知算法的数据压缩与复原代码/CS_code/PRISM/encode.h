#ifndef ENCODE_H
#define ENCODE_H

#include "Global.h"

//获取指定帧内的数据。
int GetOneFrameData(int frameIndex,FILE *fp,FrameParament &curFrameData);

//获取一帧中坐标为(x,y)的宏块数据
int GetBlockData(int x,int y,FrameParament curFrameData);

//对一帧数据进行编码。	
int EncodeIFrame(FrameParament curFrameData,int16 encodeLumaData[XX][YY]);

//对P帧数据进行编码。
int EncodePFrame(FrameParament curFrameData,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS]);

//对I帧宏块进行编码。
int EncodeMacroBlock();

//分离要进行Wyner-ziv编码的数据和进行传统编码的数据
int GetSyndromeSource();

//将P帧宏块进行编码。
int SyndromeEncode();

//设置hash校验子和伴随子等。
int SetPrismData(int x,int y,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS]);

int DCT(uint8 source[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH],int16 W[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH]);


int Quant(int x,int y, int16 W[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH]);


//将量化后的宏块保存到encodeLuma[XX][YY]中。
int SetEncodeData(int x,int y,int16 encodeLumaData[XX][YY]);

void SaveEncodeFile(int x,int y);

#endif