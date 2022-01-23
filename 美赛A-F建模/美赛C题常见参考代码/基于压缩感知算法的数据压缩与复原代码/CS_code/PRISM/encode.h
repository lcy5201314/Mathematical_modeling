#ifndef ENCODE_H
#define ENCODE_H

#include "Global.h"

//��ȡָ��֡�ڵ����ݡ�
int GetOneFrameData(int frameIndex,FILE *fp,FrameParament &curFrameData);

//��ȡһ֡������Ϊ(x,y)�ĺ������
int GetBlockData(int x,int y,FrameParament curFrameData);

//��һ֡���ݽ��б��롣	
int EncodeIFrame(FrameParament curFrameData,int16 encodeLumaData[XX][YY]);

//��P֡���ݽ��б��롣
int EncodePFrame(FrameParament curFrameData,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS]);

//��I֡�����б��롣
int EncodeMacroBlock();

//����Ҫ����Wyner-ziv��������ݺͽ��д�ͳ���������
int GetSyndromeSource();

//��P֡�����б��롣
int SyndromeEncode();

//����hashУ���ӺͰ����ӵȡ�
int SetPrismData(int x,int y,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS]);

int DCT(uint8 source[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH],int16 W[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH]);


int Quant(int x,int y, int16 W[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH]);


//��������ĺ�鱣�浽encodeLuma[XX][YY]�С�
int SetEncodeData(int x,int y,int16 encodeLumaData[XX][YY]);

void SaveEncodeFile(int x,int y);

#endif