#ifndef DECODE_H
#define DECODE_H
#include "Global.h"

//��ȡһ֡������Ϊ(i,j)�ĺ������
int GetDecodeBlockData(int x,int y,int16 encodeLumaData[XX][YY]);

//��I֡���ݽ��н��롣	
int DecodeIFrame(int16 encodeLumaData[XX][YY],FrameParament &curFrameData);

//��P֡���ݽ��б��롣
int DecodePFrame(int index,FrameParament &curFrameData,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS]);

//�Ƚ���LDPC����Ȼ��ͨ��HashУ��ó��Ƿ��Ǻ��ʵı���Ϣ��
int SyndromeDecode(int x,int y);

//P֡������QuantBlockData,encodeLumaData����;
int SetQuantBlockData(int x,int y,int16 uPrismLumaData[X][Y][UPRISMBITS],
					  int16 encodeLumaData[XX][YY]);

//��һ��������н��롣
int DecodeMacroBlock();	

int IQuantAndIDCT(int x,int y,int16 source[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH]);

//�����ĺ�鱣�浽curDecodedFrameData.LumaData[XX][YY]��;
int SetDecodeData(int x,int y,FrameParament &curFrameData);

//��һ֡�Ѿ���������ݴ浽DECODE_FILE��
int WriteOneDecodeFrame(int frameIndex,FILE *fp1,FrameParament curFrameData);

//��ȡһ�������Ӻ��
int GetSyndromeBits(int x,int y,uint8 syndromeLumaBits[X][Y][rows]);

//�õ���Ѱ�������Ϣ
int GetBestSideInfo(int x,int y,int16 syndromeCheck[X][Y*CHECKBITS]);

//LDPC����
int LDPCDecode();



//��ͨ��HashУ��ó�����Ϣ�ٽ���LDPC���롣
int SyndromeDecode3();
int SAD(int x,int y,int dx,int dy,int16 syndromeCheck[X][Y*CHECKBITS]);
int SHD(int x,int y,int dx,int dy,int16 syndromeCheck[X][Y*CHECKBITS]);
void SaveDecodeFile(int i,int j);
//LDPC����Ҫ�õĺ���
void qin2tin(double tin[rows][maxCheckDegree],double qin[cols][maxVariableDegree],double q0[cols]);
void tin2qin(double qin[cols][maxVariableDegree],double tin[rows][maxCheckDegree],uint8 *syndrome);
double Operator(double *temp, int *temp1,int j,bool flag);
void GetSyndrome(uint8 *x,uint8 *syndrome);
bool CompareSyndrome(uint8 *syndrome,uint8 *xxSyndrome);

#endif