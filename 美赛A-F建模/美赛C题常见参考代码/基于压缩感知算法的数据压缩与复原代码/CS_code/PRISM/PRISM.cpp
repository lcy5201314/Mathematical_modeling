#include "Global.h"
#include "encode.h"
#include "decode.h"

int variableNodeOnes[cols][maxVariableDegree];
int checkNodeOnes[rows][maxCheckDegree];
double theoryPbFrame[FRAME_NUM-1][X][Y];	//ͳ�Ƴ���������������;
int16 MV[FRAME_NUM-1][X][Y][2];			    //�˶�ʸ��	
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
	FrameParament curFrameData;								//���ڱ����ԭʼ�����ж�ȡ���ļ���
	int16 encodeLumaData[XX][YY];							//������û�н��зֲ�ʽ�����һ֡����
	int16 uPrismLumaData[X][Y][UPRISMBITS];					//P֡�в����а����ӱ�������					
	uint8 syndromeLumaBits[X][Y][rows];						//P֡�а����ӱ�������
	int16 syndromeCheck[X][Y*CHECKBITS];					//P֡�еĺ��У��ֵ

	FILE *fp,*fp1;											//fpָ����ƵԴ�ļ�
	if((fp = fopen(OPEN_FILE,"rb"))==NULL)	return 0;		//����ƵԴ�ļ�
	if((fp1 = fopen(DECODE_FILE,"wbt"))==NULL)	return 0;	//�򿪽����ļ�д������
	
	if((fp_enQuant = fopen(ENCODE_READ,"wt"))==NULL)	return 0;		//��������������
	if((fp_enDCT = fopen("DCT.txt","wt"))==NULL)	return 0;		//��������������
	if((fp_de = fopen(DECODE_READ,"wt"))==NULL)	return 0;		//���������ǰ����

	if((fp_mv = fopen(MV_READ,"wt"))==NULL)	return 0;		//�����˶�ʸ��
	if((fp_source = fopen("sourcex.txt","wt"))==NULL)	return 0;		//���������ǰ����
	if((fp_reSource = fopen("reSource.txt","wt"))==NULL)	return 0;		//���������ǰ����
	if((fp_sideInfo = fopen(SIDE_READ,"wt"))==NULL)	return 0;		//�����ҵ�����ѱ���Ϣ
	if ((fp_sad = fopen("sad.txt","wt"))==NULL) return 0;

	if((DecodeInfo = fopen("decodeInfo.txt","wt"))==NULL)	return 0;
	fprintf(DecodeInfo,"���ڽ������ز���\n");
	fprintf(DecodeInfo,"MAX_MOTION=%d,MAX_ITER=%d,INTRAPERIOD=%d\n",MAX_MOTION,MAX_ITER,INTRAPERIOD);

//	GetMV();
	GetFileData();

	totalEncodeTime = 0;
	totalDecodeTime = 0;
	for (i=0;i<FRAME_NUM;i++)
	{
		printf("��ʼ�����%d֡......\n",i);
		fprintf(DecodeInfo,"��ʼ�����%d֡......\n",i);
		fprintf(fp_enDCT,"��%d֡����:\n",i);
		fprintf(fp_enQuant,"��%d֡����:\n",i);

		fprintf(fp_de,"��%d֡����:\n",i);
		
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
		printf("��%d֡������ɣ�����ʱ��:%lf\n",i,encodeDuration/CLOCKS_PER_SEC);
		fprintf(DecodeInfo,"��%d֡������ɣ�����ʱ��:%lf\n",i,encodeDuration/CLOCKS_PER_SEC);

		printf("��ʼ�����%d֡......\n",i);
		fprintf(DecodeInfo,"��ʼ�����%d֡......\n",i);

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
		printf("��%d֡������ɣ�����ʱ��:%lf\n",i,decodeDuration / CLOCKS_PER_SEC);
		fprintf(DecodeInfo,"��%d֡������ɣ�����ʱ��:%lf\n",i,decodeDuration / CLOCKS_PER_SEC);

	}


	printf("����֡�������ɣ���%d֡��������ʱ��%lf,������ʱ��%lf\n",i,totalEncodeTime/CLOCKS_PER_SEC,totalDecodeTime/CLOCKS_PER_SEC);
	fprintf(DecodeInfo,"����֡�������ɣ���%d֡��������ʱ��%lf,������ʱ��%lf\n",i,totalEncodeTime/CLOCKS_PER_SEC,totalDecodeTime/CLOCKS_PER_SEC);

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
	��ȡtheoryPb.pb,variable.vc��check.vc�ļ�
		                                                                */
/************************************************************************/
bool GetFileData()
{
	//��ȡtheoryPb.pb�ļ������ͳ������������
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
	
	//��ȡvariable.vc��check.vc�ļ������variableNodeOnes,checkNodeOnes
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