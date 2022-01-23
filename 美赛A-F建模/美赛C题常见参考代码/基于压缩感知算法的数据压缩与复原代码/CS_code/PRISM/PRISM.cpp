#include "Global.h"
#include "encode.h"
#include "decode.h"

int variableNodeOnes[cols][maxVariableDegree];
int checkNodeOnes[rows][maxCheckDegree];
double theoryPbFrame[FRAME_NUM-1][X][Y];	//统计出来的理论误码率;
int16 MV[FRAME_NUM-1][X][Y][2];			    //运动矢量	
FILE *DecodeInfo;
FILE *fp_enDCT;
FILE *fp_enQuant;
FILE *fp_de;
FILE *fp_mv;
FILE *fp_source;
FILE *fp_reSource;
FILE *fp_sideInfo;
FILE *fp_sad;

bool GetFileData();
bool GetMV();

int main()
{
	int i;
	clock_t start, finish1,finish2;
	double encodeDuration,decodeDuration;
	double totalEncodeTime,totalDecodeTime;
	FrameParament curFrameData;								//用于保存从原始序列中读取的文件。
	int16 encodeLumaData[XX][YY];							//编码后的没有进行分布式编码的一帧数据
	int16 uPrismLumaData[X][Y][UPRISMBITS];					//P帧中不进行伴随子编码数据					
	uint8 syndromeLumaBits[X][Y][rows];						//P帧中伴随子编码数据
	int16 syndromeCheck[X][Y*CHECKBITS];					//P帧中的宏块校验值

	FILE *fp,*fp1;											//fp指向视频源文件
	if((fp = fopen(OPEN_FILE,"rb"))==NULL)	return 0;		//打开视频源文件
	if((fp1 = fopen(DECODE_FILE,"wbt"))==NULL)	return 0;	//打开解码文件写入数据
	
	if((fp_enQuant = fopen(ENCODE_READ,"wt"))==NULL)	return 0;		//保存量化后数据
	if((fp_enDCT = fopen("DCT.txt","wt"))==NULL)	return 0;		//保存量化后数据
	if((fp_de = fopen(DECODE_READ,"wt"))==NULL)	return 0;		//保存解量化前数据

	if((fp_mv = fopen(MV_READ,"wt"))==NULL)	return 0;		//保存运动矢量
	if((fp_source = fopen("sourcex.txt","wt"))==NULL)	return 0;		//保存解量化前数据
	if((fp_reSource = fopen("reSource.txt","wt"))==NULL)	return 0;		//保存解量化前数据
	if((fp_sideInfo = fopen(SIDE_READ,"wt"))==NULL)	return 0;		//保存找到的最佳边信息
	if ((fp_sad = fopen("sad.txt","wt"))==NULL) return 0;

	if((DecodeInfo = fopen("decodeInfo.txt","wt"))==NULL)	return 0;
	fprintf(DecodeInfo,"关于解码的相关参数\n");
	fprintf(DecodeInfo,"MAX_MOTION=%d,MAX_ITER=%d,INTRAPERIOD=%d\n",MAX_MOTION,MAX_ITER,INTRAPERIOD);

//	GetMV();
	GetFileData();

	totalEncodeTime = 0;
	totalDecodeTime = 0;
	for (i=0;i<FRAME_NUM;i++)
	{
		printf("开始编码第%d帧......\n",i);
		fprintf(DecodeInfo,"开始编码第%d帧......\n",i);
		fprintf(fp_enDCT,"第%d帧数据:\n",i);
		fprintf(fp_enQuant,"第%d帧数据:\n",i);

		fprintf(fp_de,"第%d帧数据:\n",i);
		
		start = clock();
 		GetOneFrameData(i,fp,curFrameData);
		if (i%INTRAPERIOD==0)
		{
			EncodeIFrame(curFrameData,encodeLumaData);
		}
		else
		{
			EncodePFrame(curFrameData,uPrismLumaData,syndromeLumaBits,syndromeCheck);
		}
		
		finish1 = clock();
		encodeDuration = (double)(finish1 - start);
		totalEncodeTime += encodeDuration;
		printf("第%d帧编码完成，编码时长:%lf\n",i,encodeDuration/CLOCKS_PER_SEC);
		fprintf(DecodeInfo,"第%d帧编码完成，编码时长:%lf\n",i,encodeDuration/CLOCKS_PER_SEC);

		printf("开始解码第%d帧......\n",i);
		fprintf(DecodeInfo,"开始解码第%d帧......\n",i);

		finish1 = clock();
		if (i%INTRAPERIOD==0)
		{
			DecodeIFrame(encodeLumaData,curFrameData);
		}
		else
		{
			DecodePFrame(i,curFrameData,uPrismLumaData,syndromeLumaBits,syndromeCheck);
		}

		
		WriteOneDecodeFrame(i,fp1,curFrameData);
 		finish2 = clock();
		decodeDuration = (double)(finish2 - finish1);
		totalDecodeTime +=decodeDuration;
		printf("第%d帧解码完成，解码时长:%lf\n",i,decodeDuration / CLOCKS_PER_SEC);
		fprintf(DecodeInfo,"第%d帧解码完成，解码时长:%lf\n",i,decodeDuration / CLOCKS_PER_SEC);

	}


	printf("所有帧编解码完成，总%d帧，编码总时间%lf,解码总时间%lf\n",i,totalEncodeTime/CLOCKS_PER_SEC,totalDecodeTime/CLOCKS_PER_SEC);
	fprintf(DecodeInfo,"所有帧编解码完成，总%d帧，编码总时间%lf,解码总时间%lf\n",i,totalEncodeTime/CLOCKS_PER_SEC,totalDecodeTime/CLOCKS_PER_SEC);

	fclose(fp_enDCT);
	fclose(fp_enQuant);
	fclose(fp_de);
	fclose(fp);
	fclose(fp1);
	fclose(DecodeInfo);
	fclose(fp_mv);
	fclose(fp_source);
	fclose(fp_reSource);
	fclose(fp_sideInfo);
	return 1;
}


/************************************************************************/
/* 
	读取theoryPb.pb,variable.vc和check.vc文件
		                                                                */
/************************************************************************/
bool GetFileData()
{
	//读取theoryPb.pb文件，获得统计理论误码率
	FILE *fp_pb;		 
	if((fp_pb = fopen("theoryPb_qcif.pb","rbt"))==NULL)	return false;
	fseek(fp_pb,0,SEEK_SET);
	fread(theoryPbFrame[0],sizeof(double),X*Y*(FRAME_NUM-1),fp_pb);
	fclose(fp_pb);

// 	for (int index =0;index<FRAME_NUM-1;index++)
// 	{
// 		for(int i=0;i<X;i++)
// 		{
// 			for(int j=0;j<Y;j++)
// 			{
// 				printf("theoryPbFrame[%d][%d][%d]=%lf\n",index,i,j,theoryPbFrame[index][i][j]);
// 			}
// 		}
// 		printf("\n");
// 		getch();
// 	}
	
	//读取variable.vc与check.vc文件，获得variableNodeOnes,checkNodeOnes
	FILE *fpVar,*fpChe;
	if((fpVar=fopen("variable.vc","rbt"))==NULL)	return false;
	fseek(fpVar,0,SEEK_SET);
	fread(variableNodeOnes[0],sizeof(int),cols*maxVariableDegree,fpVar);
// 	int num=fread(variableNodeOnes[0],sizeof(int),cols*maxVariableDegree,fpVar);
// 	printf("%d,%d",num,cols*maxVariableDegree);
	
	if((fpChe=fopen("check.vc","rbt"))==NULL)	return false;
	fseek(fpChe,0,SEEK_SET);
	fread(checkNodeOnes[0],sizeof(int),rows*maxCheckDegree,fpChe);
	fclose(fpVar);
	fclose(fpChe);
// 	for(int i=0;i<cols;i++)
// 	{
// 		printf("%d:",i);
// 		for(int j=0;j<maxVariableDegree;j++)
// 		{
// 			printf("%d ",variableNodeOnes[i][j]);
// 		}
// 		printf("\n");
// 	}
// 	getch();
// 
// 	for(i=0;i<rows;i++)
// 	{
// 		printf("%d:",i);
// 		for(int j=0;j<maxCheckDegree;j++)
// 		{
// 			printf("%d ",checkNodeOnes[i][j]);
// 
// 		}
// 		printf("\n");
// 	}
	return true;
}


bool GetMV()
{
	FILE *fp_mv;
	if((fp_mv = fopen("Motion_qcif.mv","rbt"))==NULL) return false;

	for (int index = 0;index<FRAME_NUM-1;index++)
	{
		fseek(fp_mv,X*Y*2*sizeof(int16)*index,SEEK_SET);
		fread(MV[index],sizeof(int16),X*Y*2,fp_mv);
	}
	fclose(fp_mv);

// 	for (index = 0;index<FRAME_NUM-1;index++)
// 	{
// 		printf("index=%d\n",index);
// 		for( int i=0;i<X;i++)
// 		{
// 			for( int j=0;j<Y;j++)
// 			{
// 				printf("(%d,%d) ",MV[index][i][j][0],MV[index][i][j][1]);
// 			}
// 			printf("\n");
// 		}
// 		getch();
// 	}
// 
	return true;
}