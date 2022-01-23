#include "decode.h"


uint32 residual;
bool finish = false;
uint8 decodeBlockBits[cols];		//���ڼ�¼�����ӽ�����ĺ������
int16 sideInfoSource[PRISMBITS];	//����Ӧ�ı���Ϣʮ���Ʊ�ʾ
uint8 sideInfo[cols];				//����Ӧ�ı���Ϣ�����Ʊ�ʾ
int16 sideLumaData[XX][YY];			//���ڴ�Ž����һ֡������������һ֡�����ӽ���ı���Ϣ
double theoryPb = 0;				//���ֵ��Ҫͳ�Ʋ���ȷ����
uint32 sad = 0;

/************************************************************************/
/* 
   ��I֡���ݽ��н��롣
                                                                        */
/************************************************************************/
int DecodeIFrame(int16 encodeLumaData[XX][YY],FrameParament &curFrameData)
{
	int i,j;
	
	//��encodeLumaData����д��sideLumaData��,��Ϊ��һ֡����ı���Ϣ
	memcpy(sideLumaData,encodeLumaData,XX*YY*sizeof(int16));

	for(i=0;i<X;i++)
	{
		for(j=0;j<Y;j++)
		{
			GetDecodeBlockData(i, j, encodeLumaData);
			SaveDecodeFile(i,j);
			DecodeMacroBlock();
			SetDecodeData(i,j,curFrameData); //�������ĺ�鱣�浽curFrameData.LumaData[XX][YY]
		}
	}
	return 0;
}

/************************************************************************/
/* 
   ��P֡���ݽ��н��롣
                                                                        */
/************************************************************************/
int DecodePFrame(int index,FrameParament &curFrameData,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS])
{
	int16 encodeLumaData[XX][YY] = {0};//���ڴ�ŵ�ǰ֡δ���з�������DCT���������ݡ�
	int i,j;

	for(i=0;i<X;i++)
	{
		for(j=0;j<Y;j++)
		{
			GetSyndromeBits(i,j,syndromeLumaBits);
			GetBestSideInfo(i,j,syndromeCheck);//��ͨ��HashУ��ó�����Ϣ�ٽ���LDPC���롣
			theoryPb = 0.362296;// theoryPbFrame[index][i][j]; //0.183233;//
			SyndromeDecode3();
			SetQuantBlockData(i,j,uPrismLumaData,encodeLumaData);
			SaveDecodeFile(i,j);
			DecodeMacroBlock();
			SetDecodeData(i,j,curFrameData); //�������ĺ�鱣�浽curFrameData.LumaData[XX][YY]
		}
		fprintf(fp_mv,"\n");
	}

	//����ǰ֡�������LumaData����д��sideLumaData��,��Ϊ��һ֡����ı���Ϣ
	memcpy(sideLumaData,encodeLumaData,XX*YY*sizeof(int16));

	return 0;
}

/************************************************************************/
/*
	��ȡһ����01���б�ʾ�İ����Ӻ�����С�
					                                                    */
/************************************************************************/
int GetSyndromeBits(int x,int y,uint8 syndromeLumaBits[X][Y][rows])
{	
	for(int k=0;k<rows;k++)
	{
		syndromeBlockBits[k] = syndromeLumaBits[x][y][k];
//		printf("%d ",syndromeLumaBits[x][y][k]);
	}

// 	getch();
	return 1;
}

/************************************************************************/
/*
	�õ���Ѱ�������Ϣ
					                                                    */
/************************************************************************/
int GetBestSideInfo(int x,int y,int16 syndromeCheck[X][Y*CHECKBITS])
{	

	int i,j;
	int dx,dy;
	int k=0;
	sad=0xffffff;
	residual = 1000;
	finish = false;
	dx = 0;
	dy = 0;
 
// 	static bool writesidefile = true;
// 	if (writesidefile)
// 	{
// 		for (i=0;i<X;i++)
// 		{
// 			for (j=0;j<Y;j++)
// 			{
// 				int ox = i<<4;
// 				int oy = j<<4;
// 				fprintf(fp_sideInfo,"i=%d,j=%d\n",i,j);
// 				for (int k=0;k<BLOCK_HEIGHT;k++)
// 				{
// 					for (int m=0;m<BLOCK_WIDTH;m++)
// 					{
// 						fprintf(fp_sideInfo,"%d ",sideLumaData[ox+k][oy+m]);
// 					}
// 					fprintf(fp_sideInfo,"\n");
// 				}	
// 
// 			}
// 		}
// 		writesidefile = false;
// 	}
 	

	for(i=-MAX_MOTION;i<=MAX_MOTION;i++)
	{
		for(j=-MAX_MOTION;j<=MAX_MOTION;j++)
		{
			if (SAD(x,y,i,j,syndromeCheck)==1) //ͨ��SAD��ѡ�����Ϣ
			//if (SHD(x,y,i,j,syndromeCheck)==1)	 //ͨ����������ѡ�����Ϣ
			{
				dx = i;
				dy = j;
				fprintf(fp_sad,"i=%d,j=%d,dx=%d,dy=%d,sad=%d\n",x,y,dx,dy,sad);
			}
		}
	}

	fprintf(fp_mv,"(%d,%d) ",dx,dy);

	int rx,ry,ox,oy;
	ox = x<<4;			// x*BLOCK_HEIGHT
	oy = y<<4;			// y*BLOCK_WIDTH	
	rx = dx + ox;
	ry = dy + oy;
	
	//�����ĺ�����򳬳�һ֡��С�ı߽硣
	if( rx<0 || ry<0 || (rx+BLOCK_HEIGHT)>XX || (ry+BLOCK_WIDTH)>YY)
	{
//		fprintf(fp_sideInfo,"i=%d,j=%d,error!\n",x,y);
		return 0;
	}
	//�õ���(rx,ry)Ϊ�󶥵�����ı���Ϣ����ʮ���Ʊ�ʾ
	int count = 0;
	int tempIndex = 0;
	int synsourceIndex = 0;
	int subox = 0;
	int suboy = 0;

	for (i=0;i<subX;i++)
	{
		subox = i<<2;    //i*SUBBLOCK_HEIGHT;
		for (j=0;j<subY;j++)
		{

			suboy = j<<2;    //j*SUBBLOCK_WIDTH;	
			synsourceIndex = tempIndex<<2;
			sideInfoSource[synsourceIndex+0] = sideLumaData[rx+subox+0][ry+suboy+0];
			sideInfoSource[synsourceIndex+1] = sideLumaData[rx+subox+0][ry+suboy+1];
			sideInfoSource[synsourceIndex+2] = sideLumaData[rx+subox+1][ry+suboy+0];
			sideInfoSource[synsourceIndex+3] = sideLumaData[rx+subox+2][ry+suboy+0];
			tempIndex++;
		}
	}
	
	fprintf(fp_sideInfo,"i=%d,j=%d\n",x,y);
	tempIndex = 0;
	for (i=0;i<subX;i++)
	{
		for (j=0;j<subY;j++)
		{
			synsourceIndex = tempIndex<<2;
			for (k=0;k<4;k++)
			{
				fprintf(fp_sideInfo,"%4d ",sideInfoSource[synsourceIndex+k]);
			}
			tempIndex++;
		}
	}
	fprintf(fp_sideInfo,"\n");


	//��sideInfoSource�����������ת����16bits�Ķ����������浽����sideInfo��
	tempIndex = 0;
	count = 0;
	for (i=0;i<PRISMBITS;i++)
	{
//		printf("%d\n",sideInfoSource[i]);
		tempIndex = count<<4;		//ͨ������������������䣬���ټ��㡣
		count++;
		for (k=0;k<16;k++)
		{
			if (sideInfoSource[i] & 1<<(15-k))
			{
				sideInfo[tempIndex+k] = 1;
			}
			else
			{
				sideInfo[tempIndex+k] = 0;
			}
//			printf("%d ",sideInfo[tempIndex+k]);

		}
//		printf("\n");
//		getch();
	}

// 	sideInfo[cols-1] = 0;
// 	sideInfo[cols-2] = 0;

//	getch();
//	printf("//////////////////////////////////////////////////////////////////////\n");
	return 1;
}

/************************************************************************/
/* 
	��һ֡�����еĺ����а����ӽ���,����Ϣ��sideLumaData[XX][YY]�����С�
	�����ӽ��������ݴ����syndromeSource���顣
	                                                                    */
/************************************************************************/
int SyndromeDecode3()
{
	bool test =false;
	int count = 0;
	int tempIndex = 0;

	LDPCDecode();

	/*����LDPC���뻹ԭ�Ķ���������*/
	for (int k=0;k<cols;k++)
	{
		fprintf(fp_reSource,"%d ",decodeBlockBits[k]);
	}
	fprintf(fp_reSource,"\n");

	tempIndex = 0;
	count = 0;
	for (int i=0;i<PRISMBITS;i++)
	{
		tempIndex = count<<4;		//ͨ������������������䣬���ټ��㡣
		count++;
		syndromeSource[i] = 0;
		for (int j=0;j<16;j++)
		{
			syndromeSource[i] += ((decodeBlockBits[tempIndex+(15-j)])<<j);
		}
	}



	return 1;
}

/************************************************************************/
/* 
	P֡������QuantBlockData��encodeLumaData����
	                                                                    */
/************************************************************************/
int SetQuantBlockData(int x,int y,int16 uPrismLumaData[X][Y][UPRISMBITS],int16 encodeLumaData[XX][YY])
{
	int tempIndex = 0;
	int ox = x<<4;
	int oy = y<<4;
	int subox = 0;
	int suboy = 0; 
	int synsourceIndex = 0;
	int blockSourceIndex = 0;

	for (int i=0;i<subX;i++)
	{
		subox = i<<2;    //x*SUBBLOCK_HEIGHT;
		for (int j=0;j<subY;j++)
		{
			suboy = j<<2;    //y*SUBBLOCK_WIDTH;		
			synsourceIndex = tempIndex<<2;
			blockSourceIndex = tempIndex*12 - 4;
			QuantBlockData[subox+0][suboy+0] = syndromeSource[synsourceIndex+0];
			encodeLumaData[ox+subox+0][oy+suboy+0] = QuantBlockData[subox+0][suboy+0];
			
			QuantBlockData[subox+0][suboy+1] = syndromeSource[synsourceIndex+1];
			encodeLumaData[ox+subox+0][oy+suboy+1] = QuantBlockData[subox+0][suboy+1];

			QuantBlockData[subox+1][suboy+0] = syndromeSource[synsourceIndex+2];
			encodeLumaData[ox+subox+1][oy+suboy+0] = QuantBlockData[subox+1][suboy+0];

			QuantBlockData[subox+2][suboy+0] = syndromeSource[synsourceIndex+3];
			encodeLumaData[ox+subox+2][oy+suboy+0] = QuantBlockData[subox+2][suboy+0];

			QuantBlockData[subox+1][suboy+1] = uPrismLumaData[x][y][blockSourceIndex+4];
			encodeLumaData[ox+subox+1][oy+suboy+1] = QuantBlockData[subox+1][suboy+1];

			QuantBlockData[subox+0][suboy+2] = uPrismLumaData[x][y][blockSourceIndex+5];
			encodeLumaData[ox+subox+0][oy+suboy+2] = QuantBlockData[subox+0][suboy+2];
			
			QuantBlockData[subox+0][suboy+3] = uPrismLumaData[x][y][blockSourceIndex+6];
			encodeLumaData[ox+subox+0][oy+suboy+3] = QuantBlockData[subox+0][suboy+3];
			
			QuantBlockData[subox+1][suboy+2] = uPrismLumaData[x][y][blockSourceIndex+7];
			encodeLumaData[ox+subox+1][oy+suboy+2] = QuantBlockData[subox+1][suboy+2];
			
			QuantBlockData[subox+2][suboy+1] = uPrismLumaData[x][y][blockSourceIndex+8];
			encodeLumaData[ox+subox+2][oy+suboy+1] = QuantBlockData[subox+2][suboy+1];
			
			QuantBlockData[subox+3][suboy+0] = uPrismLumaData[x][y][blockSourceIndex+9];
			encodeLumaData[ox+subox+3][oy+suboy+0] = QuantBlockData[subox+3][suboy+0];
			
			QuantBlockData[subox+3][suboy+1] = uPrismLumaData[x][y][blockSourceIndex+10];
			encodeLumaData[ox+subox+3][oy+suboy+1] = QuantBlockData[subox+3][suboy+1];
			
			QuantBlockData[subox+2][suboy+2] = uPrismLumaData[x][y][blockSourceIndex+11];
			encodeLumaData[ox+subox+2][oy+suboy+2] = QuantBlockData[subox+2][suboy+2];
			
			QuantBlockData[subox+1][suboy+3] = uPrismLumaData[x][y][blockSourceIndex+12];
			encodeLumaData[ox+subox+1][oy+suboy+3] = QuantBlockData[subox+1][suboy+3];
			
			QuantBlockData[subox+2][suboy+3] = uPrismLumaData[x][y][blockSourceIndex+13];
			encodeLumaData[ox+subox+2][oy+suboy+3] = QuantBlockData[subox+2][suboy+3];
			
			QuantBlockData[subox+3][suboy+2] = uPrismLumaData[x][y][blockSourceIndex+14];
			encodeLumaData[ox+subox+3][oy+suboy+2] = QuantBlockData[subox+3][suboy+2];
			
			QuantBlockData[subox+3][suboy+3] = uPrismLumaData[x][y][blockSourceIndex+15];
			encodeLumaData[ox+subox+3][oy+suboy+3] = QuantBlockData[subox+3][suboy+3];
		
			tempIndex++;
		}
	}

	
	return 1;
}

/************************************************************************/
/* 
	LDPC���벿�֡�
	                                                                    */
/************************************************************************/
int LDPCDecode()
{
	int k=0;
	//ͨ������������syndromeBlockBits�����ϢsideInfo���н��롣
	int iter;
	bool sign = false;
	uint8 xxSyndrome[rows] = {0};
	double q0[cols] = {0};
	double qin[cols][maxVariableDegree] = {0};
	double tin[rows][maxCheckDegree] = {0};
	
                   
	/*��ʼ��q0*/
	double tempq0=log((1-theoryPb)/theoryPb)/log(2);
	for(k=0;k<cols;k++)
	{
		q0[k]=(1-2*sideInfo[k])*tempq0;
	}
	
	/*ѭ�����룻*/
	for(iter=0;iter<MAX_ITER;iter++)
	{
//		printf("i=%d,j=%d,x=%d,y=%d,iter=%d,residual=%d,tempResidual=%d\n",i,j,x,y,iter,residual,tempRisual);

// 		fprintf(DecodeInfo,"i=%d,j=%d,x=%d,y=%d,iter=%d,residual=%d\n",i,j,x,y,iter,residual);
		qin2tin(tin,qin,q0);
		tin2qin(qin,tin,syndromeBlockBits);
	
		/*�о�;*/
		for(int l=0;l<cols;l++)
		{
			double temp=0;
			double tempXX=0;
			for(k=0;k<maxVariableDegree;k++)
			{
				if(variableNodeOnes[l][k]==-1) 
				{
					continue;
				}

				temp+=qin[l][k];
			}

			tempXX=q0[l]+temp;
			
			if(tempXX>=0)
			{
				decodeBlockBits[l]=0;
			}
			else
			{
				decodeBlockBits[l]=1;
			}
		}
		
		GetSyndrome(decodeBlockBits,xxSyndrome);
		sign=CompareSyndrome(syndromeBlockBits,xxSyndrome);
	
		if(sign)   /*sign���ڱ��x��xx�Ƿ����,ͨ�����������ж�;true��ʾ��ȫ��ȷ����;*/
		{
//			printf("iter=%d,��ȷ����\n",iter);
			finish = true;
			return 1;
		}
	}
	
	return 1;
}




/************************************************************************/
/*
	��ȡһ���Ѿ�����ĺ��׼������
					                                                    */
/************************************************************************/
int GetDecodeBlockData(int x,int y,int16 encodeLumaData[XX][YY])
{	
	int ox = x<<4;                //i*BLOCK_HEIGHT;
	int oy = y<<4;			      //j*BLOCK_WIDTH; 
	for (int i=0;i<BLOCK_HEIGHT;i++)
	{
		for (int j=0;j<BLOCK_WIDTH;j++)
		{
			QuantBlockData[i][j] = encodeLumaData[ox+i][oy+j];
		}
	}
	return 1;
}


/************************************************************************/
/*
	��һ��������н��룬�����뷴�����ͷ�DCT����浽MacroBlockData��
					                                                    */
/************************************************************************/
int DecodeMacroBlock()
{
	int i,j;
	int16 subMacroBLock[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH];
	for (i=0;i<subX;i++)
	{
		for (j=0;j<subY;j++)
		{
			for (int k=0;k<SUBBLOCK_HEIGHT;k++)
			{
				for (int l=0;l<SUBBLOCK_WIDTH;l++)
				{
					subMacroBLock[k][l] = QuantBlockData[(i<<2)+k][(j<<2)+l];
				}
			}
			IQuantAndIDCT(i,j,subMacroBLock);
		}
	}

	return 1;
}



/************************************************************************/
/* 
	��4x4�ӿ���з������뷴����DCT�任
	                                                                    */
/************************************************************************/
int IQuantAndIDCT(int x,int y,int16 source[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH])
{
	int i,j;
	int ox = x<<2;         // x*SUBBLOCK_HEIGHT;
	int oy = y<<2;         // y*SUBBLOCK_WIDTH;
	int16 temp[2][2]={0};
	int qp_per,qp_rem;

	int16 xx[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH] = {0};
	int16 W[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH] = {0};

	qp_per    = (QP-MIN_QP)/6;	//��6������������
	qp_rem    = (QP-MIN_QP)%6;  //��6������������

	
	//����������
	for (i=0;i<SUBBLOCK_HEIGHT;i++)
	{
		for (j=0;j<SUBBLOCK_WIDTH;j++)
		{
			W[i][j] = source[i][j]*dequant_coef[qp_rem][i][j]*(1<<qp_per);
		}
	}
	
	//����άDCT��ת��Ϊ����һά�任�Ƚ����б任
    for (i=0;i<SUBBLOCK_WIDTH;i++)
    {
		temp[0][0] = W[0][i] + W[2][i];
		temp[0][1] = W[0][i] - W[2][i];
		temp[1][0] = W[1][i] + (W[3][i]>>1);
		temp[1][1] = (W[1][i]>>1) - W[3][i];

		xx[0][i] = temp[0][0] + temp[1][0];
		xx[1][i] = temp[0][1] + temp[1][1];
		xx[2][i] = temp[0][1] - temp[1][1];
		xx[3][i] = temp[0][0] - temp[1][0];
    }
	
	//�б任
	for (i=0;i<SUBBLOCK_HEIGHT;i++)
    {
		temp[0][0] = xx[i][0] + xx[i][2];
		temp[0][1] = xx[i][0] - xx[i][2];
		temp[1][0] = xx[i][1] + (xx[i][3]>>1);
		temp[1][1] = (xx[i][1]>>1) - xx[i][3];

		//�˴�Ӧ���ٳ�64���ο���Ƥ��Transform&�������ֵ���ڽ��벿�֡�
		MacroBlockData[i+ox][0+oy] = (temp[0][0] + temp[1][0])>>6;   
		MacroBlockData[i+ox][1+oy] = (temp[0][1] + temp[1][1])>>6;
		MacroBlockData[i+ox][2+oy] = (temp[0][1] - temp[1][1])>>6;
		MacroBlockData[i+ox][3+oy] = (temp[0][0] - temp[1][0])>>6;
   }
	return 1;
}


/************************************************************************/
/* 
	�����ĺ�鱣�浽curFrameData.LumaData[XX][YY]��;
	                                                                    */
/************************************************************************/
int SetDecodeData(int x,int y,FrameParament &curFrameData)
{
	int ox = x<<4;    //x*BLOCK_HEIGHT;
	int oy = y<<4;    //y*BLOCK_WIDTH;
	int i,j;
	for (i=0;i<BLOCK_HEIGHT;i++)
	{
		for (j=0;j<BLOCK_WIDTH;j++)
		{
			curFrameData.LumaData[i+ox][j+oy] = MacroBlockData[i][j];
		}
	}
	return 1;
}


/************************************************************************/
/* 
	��һ֡�Ѿ���������ݴ浽DECODE_FILE�ļ���
                                                                        */
/************************************************************************/
int WriteOneDecodeFrame(int frameIndex,FILE *fp1,FrameParament curFrameData) 
{
	//�洢ENCODE_FILE�ļ�
	fseek(fp1,(XX*YY*3*frameIndex)>>1,SEEK_SET);
	if(fwrite(curFrameData.LumaData[0],sizeof(uint8),XX*YY,fp1)==0) return 0;
	if(fwrite(curFrameData.UData[0],sizeof(uint8),(XX*YY)>>2,fp1)==0) return 0;
	if(fwrite(curFrameData.VData[0],sizeof(uint8),(XX*YY)>>2,fp1)==0) return 0;
	return 1;
}


/************************************************************************/
/* 
	LDPC����Ҫ�õĺ���
                                                                        */
/************************************************************************/
void qin2tin(double tin[rows][maxCheckDegree],double qin[cols][maxVariableDegree],double q0[cols])
{
	int i,j,k,l,m;
	double sum=0;
	double **qout;
/*	double qout[cols][maxVariableDegree]={0};*/
	bool flag=true;             /*���ñ�־λ����Operator�Ǽӻ��ǳ�;true��Ӧ"��"��false��Ӧ�ˣ�*/
	double *temp;
	int *tempVariable;

	/************************************************************************/
	/* Ϊqout����ռ�                                                       */
	qout=(double **)malloc(cols*sizeof(double *));
	for (i=0;i<cols;i++)
	{
		qout[i]=(double *)malloc(maxVariableDegree*sizeof(double));
	}
	for (i=0;i<cols;i++)
	{
		for (j=0;j<maxVariableDegree;j++)
		{
			qout[i][j]=0;
		}
	}

	/************************************************************************/

	for(i=0;i<cols;i++)    /*��qin�õ�qout*/
	{
		temp=qin[i];
		tempVariable=variableNodeOnes[i];
		for(j=0;j<maxVariableDegree;j++)
		{
			if(variableNodeOnes[i][j]==-1)
			{
				break;
			}

			sum=Operator(temp,tempVariable,j,flag);
			qout[i][j]=q0[i]+sum;
		}
		temp=NULL;
	}
	
/*	free(temp);*/
	for(k=0;k<rows;k++)		/*��qout�õ�tin	*/
	{
		for(l=0;l<maxCheckDegree;l++)
		{
			if(checkNodeOnes[k][l]==-1)
			{
				break;
			}
			int variableNodeIndex=checkNodeOnes[k][l];
			for(m=0;m<maxVariableDegree;m++)   
			{
				/*�õ�tin[i][j]�е�i��qout�ж�Ӧ����ţ�*/
				if(variableNodeOnes[variableNodeIndex][m]==k)
				{
					break;
				}
			}

			tin[k][l]=qout[variableNodeIndex][m];
		}

	}

	/************************************************************************/
	/* �ͷŶ�̬������ڴ�                                                   */
	for (i=0;i<cols;i++)
	{
		free(qout[i]);
	}
	free(qout);
	/************************************************************************/

}

void tin2qin(double qin[cols][maxVariableDegree],double tin[rows][maxCheckDegree],uint8 *syndrome)
{
	int i,j,k,l,m;
	double **tout;
	double prod=1;
	double tempTout;
	bool flag;
	double *temp;
	int *tempCheck;
	
	flag=false;
	/************************************************************************/
	/* ��ʼ��tout                                                           */
	tout=(double **)malloc(rows*sizeof(double *));
	for (i=0;i<rows;i++)
	{
		tout[i]=(double *)malloc(maxCheckDegree*sizeof(double));
	}

	for (i=0;i<rows;i++)
	{
		for (j=0;j<maxCheckDegree;j++)
		{
			tout[i][j]=0;
		}
	}
	/************************************************************************/
	
	
	for(i=0;i<rows;i++)
	{
		temp=tin[i];
		tempCheck=checkNodeOnes[i];
		for(j=0;j<maxCheckDegree;j++)
		{
			if(checkNodeOnes[i][j]==-1)
			{
				break;
			}

			prod=Operator(temp,tempCheck,j,flag);
			tempTout=(1-2*syndrome[i])*prod;
/*			tempTout=(1-2*syndrome[i])*pow(10,prod);*/
			tout[i][j]=log((1+tempTout)/(1-tempTout));		/*atanh(x)=(1/2)*ln((1+z)/(1-z));*/
		}
		temp=NULL;
	}

/*	free(temp);*/
	for(k=0;k<cols;k++)
	{
		for(l=0;l<maxVariableDegree;l++)
		{
			if(variableNodeOnes[k][l]==-1)
			{
				break;
			}
			int checkNodeIndex=variableNodeOnes[k][l];
			for(m=0;m<maxCheckDegree;m++)
			{
				if(checkNodeOnes[checkNodeIndex][m]==k)
				{
					break;
				}
			}

			qin[k][l]=tout[checkNodeIndex][m];
		}
	}
	/************************************************************************/
	/* �ͷŶ�̬������ڴ�                                                   */
	for (i=0;i<rows;i++)
	{
		free(tout[i]);
	}
	free(tout);
	/************************************************************************/
}


/************************************************************************/
/* ����flag��ѡ����Ӧ������                                             */
/************************************************************************/
double Operator(double *temp, int *temp1,int j,bool flag)
{
	int k;
	double tempSum;
	
	tempSum=0;

	if(flag)		/*�����Ϊ��*/
	{
		for(k=0;k<maxVariableDegree;k++)
		{
			if(k==j || *(temp1+k)==-1) continue;
			tempSum+=*(temp+k);
		}
	}
	else
	{
		tempSum=1;

		for(k=0;k<maxCheckDegree;k++)
		{
			if(k==j || *(temp1+k)==-1) continue;
			tempSum*=tanh((*(temp+k))/2);

/*�����Ϊ��;ͨ��ȡ�����������Ϊ��:ln(a*b)=ln(a)+ln(b)
			tempSum+=log(tanh((*(temp+k))/2));	
*/		}
		
	}
	return tempSum;
}


/***************************************************************************/
/* CompareSyndrome()���ڱȽ���������֮�Ƿ���ȣ���ȷ���true�����򷵻�false*/
/***************************************************************************/
bool CompareSyndrome(uint8 *syndrome,uint8 *xxSyndrome)
{
	int i=0;

	for(i=0;i<rows;i++)
	{
		if(syndrome[i]!=xxSyndrome[i])
		{
// 			printf("(%d,%d) ",syndrome[i],xxSyndrome[i]);
// 			getch();
			return false;
		}
	}
	return true;
}

/************************************************************************/
/* GetSyndrome()���������,��������syndromeָ����ָ������               */
/************************************************************************/
void GetSyndrome(uint8 *x,uint8 *syndrome)
{
	int i,j;
	int tempSyndrome=0;

	for(i=0;i<rows;i++)
	{
		tempSyndrome=0;
		for(j=0;j<maxCheckDegree;j++)
		{
			if(checkNodeOnes[i][j]==-1)
			{
				break;
			}
			tempSyndrome+=x[checkNodeOnes[i][j]];
		}
		         
		syndrome[i]=tempSyndrome%2;
	}
}




void SaveDecodeFile(int i,int j)
{
	fprintf(fp_de,"i=%d,j=%d\n",i,j);
	for (i=0;i<BLOCK_HEIGHT;i++)
	{
		for (j=0;j<BLOCK_WIDTH;j++)
		{
			fprintf(fp_de,"%4d ",QuantBlockData[i][j]);
		}
		fprintf(fp_de,"\n");
	}

}




int SAD(int x,int y,int dx,int dy,int16 syndromeCheck[X][Y*CHECKBITS])
{
	int ox=x*BLOCK_HEIGHT,oy=y*BLOCK_WIDTH;
	int rx=ox+dx,ry=oy+dy;
	if ( abs(dx)>MAX_MOTION || abs(dy)>MAX_MOTION)
		return 0;
	if( rx<0 || ry<0 || rx+BLOCK_HEIGHT>XX || ry+BLOCK_WIDTH>YY )
		return 0;
	uint32 tempSad=0;
	for (int k=0;k<CHECKBITS;k++)
	{
		int m,n;
		m = k>>4;		//k��BLOCK_WIDTH����������
		n = k - (m<<4); //k��BLOCK_WIDTH����������
		tempSad+=abs(syndromeCheck[x][y*CHECKBITS+k] - sideLumaData[rx+m][ry+n]);
	}

	if(tempSad >= sad)
		return 0;
	sad = tempSad;
	return 1;
}


int checkHamDistance(int16 a,int16 b)
{
	int temp = 0;
	for (int i = 0;i<16;i++)
	{
		if ((a & (1<<i))!=(b & (1<<i)))
		{
			temp++;
		}
	}
	return temp;

}

int SHD(int x,int y,int dx,int dy,int16 syndromeCheck[X][Y*CHECKBITS])
{
	int ox=x*BLOCK_HEIGHT,oy=y*BLOCK_WIDTH;
	int rx=ox+dx,ry=oy+dy;
	if ( abs(dx)>MAX_MOTION || abs(dy)>MAX_MOTION)
		return 0;
	if( rx<0 || ry<0 || rx+BLOCK_HEIGHT>XX || ry+BLOCK_WIDTH>YY )
		return 0;
	uint32 tempSad=0;
	for (int k=0;k<CHECKBITS;k++)
	{
		int m,n;
		m = k>>4;		//k��BLOCK_WIDTH����������
		n = k - (m<<4); //k��BLOCK_WIDTH����������

		tempSad+=checkHamDistance(syndromeCheck[x][y*CHECKBITS+k],sideLumaData[rx+m][ry+n]);
	
	}
	if(tempSad >= sad)
		return 0;
	sad = tempSad;
	return 1;
}