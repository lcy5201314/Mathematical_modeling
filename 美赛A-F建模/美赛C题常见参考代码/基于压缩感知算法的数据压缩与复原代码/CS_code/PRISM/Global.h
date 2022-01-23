#ifndef GLOBAL_H

#define GLOBAL_H

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <conio.h>
#include <memory.h>
#include <stdlib.h>

typedef int  int32;
typedef short  int16;
typedef char  int8;
typedef unsigned short  uint16;
typedef unsigned int  uint32;
typedef unsigned char  uint8;


#define CIF 1
#define QCIF 0

#if QCIF
#define XX 288
#define YY 352
#define FRAME_NUM 3
#define OPEN_FILE "mobile_cif.yuv"		   //300֡
#define DECODE_FILE "DecodeData_cif.yuv"
#define ENCODE_FILE "EncodeData_cif.zx"
#define ENCODE_READ "EncodeDataForRead_cif.txt"
#define DECODE_READ "DecodeDataForRead_qcif.txt"
#else
#define XX 144
#define YY 176
#define FRAME_NUM 3
#define OPEN_FILE "foreman_qcif.yuv"   //3֡
//#define OPEN_FILE "salesman_qcif.yuv"   //100֡
#define DECODE_FILE "DecodeData_qcif.yuv"
#define ENCODE_FILE "EncodeData_qcif.zx"
#define ENCODE_READ "EncodeDataForRead_qcif.txt"
#define DECODE_READ "DecodeDataForRead_qcif.txt"
#define MV_READ		"mv_qcif.txt"
#define SIDE_READ	"SideIndo.txt"
#endif

// #define EncodeDEBUG 0
// #define DecodeDEBUG 0

#define BLOCK_HEIGHT 16
#define BLOCK_WIDTH 16
#define SUBBLOCK_HEIGHT 4
#define SUBBLOCK_WIDTH 4

#define MAX_MOTION 16
#define MAX_ITER   20
#define MINIERROR  10	//LDPC����ѡ�����Ϣʱ����С������

#define INTRAPERIOD 2	//I֡����

#define PRISMNUMS   4   //ÿ4x4��С�����Wyner-Ziv����ĸ���
#define PRISMBITS   64	//ÿ�������а����ӱ�������ݣ�ÿ��4x4С��Ϊ4��
#define CHECKBITS	32	//����CHECKBITS��У��λ����Ϊһ������hashֵ
#define UPRISMBITS	192 //ÿ������в����а����ӱ�������ݣ�ÿ��4x4С��Ϊ12��
#define QP 5			//��������(Quanzization parameter)
#define MIN_QP 0
#define MAX_QP 51
#define Q_BITS 15
#define Qstep 1.125     //������������QPͨ�����õ�Y

#define X (XX/BLOCK_HEIGHT)      //����֡�ڸ߶�����X�����(���Ⱥ�ɫ�ȶ�һ��)
#define Y (YY/BLOCK_WIDTH)		//����֡�ڿ������Y�����
#define subX (BLOCK_HEIGHT/SUBBLOCK_HEIGHT)
#define subY (BLOCK_WIDTH/SUBBLOCK_WIDTH)

//���������Ҫ����LDPC����룻
// #define maxVariableDegree 3
// #define maxCheckDegree 6
// #define rows (PRISMBITS<<3)
// #define cols (PRISMBITS<<4)


#define maxVariableDegree 7
#define maxCheckDegree 7
#define rows 512
#define cols 1024

extern int variableNodeOnes[cols][maxVariableDegree];
extern int checkNodeOnes[rows][maxCheckDegree];
extern double theoryPbFrame[FRAME_NUM-1][X][Y];	//ͳ�Ƴ���������������;
extern int16 MV[FRAME_NUM-1][X][Y][2];			//�˶�ʸ��

typedef struct
{
	uint8 LumaData[XX][YY],UData[XX>>1][YY>>1],VData[XX>>1][YY>>1];
}FrameParament;


//����һЩȫ�ֱ���
extern uint8 MacroBlockData[BLOCK_HEIGHT][BLOCK_WIDTH];	//ԭʼ�����еĺ��
extern int16 QuantBlockData[BLOCK_HEIGHT][BLOCK_WIDTH];	//ԭʼ���к�龭��DCT������
extern int16 syndromeSource[PRISMBITS];					//Ҫ���а����ӱ���ĺ������
extern int16 BlockSource[UPRISMBITS];					//Ҫ���д�ͳ����ĺ������

extern uint8 syndromeBlockBits[rows];					//���İ��������С�

extern FILE *DecodeInfo;
extern FILE *fp_enDCT;
extern FILE *fp_enQuant;
extern FILE *fp_de;
extern FILE *fp_mv;
extern FILE *fp_source;
extern FILE *fp_reSource;
extern FILE *fp_sideInfo;
extern FILE *fp_sad;

//�����ͷ�����������ԴH.264Դ��JM8.6��
static const int quant_coef[6][4][4] = {
  {{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243},{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243}},
  {{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660},{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660}},
  {{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194},{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194}},
  {{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647},{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647}},
  {{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355},{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355}},
  {{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893},{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893}}
};

static const int dequant_coef[6][4][4] = {
  {{10, 13, 10, 13},{ 13, 16, 13, 16},{10, 13, 10, 13},{ 13, 16, 13, 16}},
  {{11, 14, 11, 14},{ 14, 18, 14, 18},{11, 14, 11, 14},{ 14, 18, 14, 18}},
  {{13, 16, 13, 16},{ 16, 20, 16, 20},{13, 16, 13, 16},{ 16, 20, 16, 20}},
  {{14, 18, 14, 18},{ 18, 23, 18, 23},{14, 18, 14, 18},{ 18, 23, 18, 23}},
  {{16, 20, 16, 20},{ 20, 25, 20, 25},{16, 20, 16, 20},{ 20, 25, 20, 25}},
  {{18, 23, 18, 23},{ 23, 29, 23, 29},{18, 23, 18, 23},{ 23, 29, 23, 29}}
};

#endif