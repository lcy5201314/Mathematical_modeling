#include "decode.h"


uint32 residual;
bool finish = false;
uint8 decodeBlockBits[cols];		//用于记录伴随子解码出的宏块数据
int16 sideInfoSource[PRISMBITS];	//宏块对应的边信息十进制表示
uint8 sideInfo[cols];				//宏块对应的边信息二进制表示
int16 sideLumaData[XX][YY];			//用于存放解码后一帧的数据用于下一帧伴随子解码的边信息
double theoryPb = 0;				//这个值需要统计才能确定。
uint32 sad = 0;

/************************************************************************/
/* 
   对I帧数据进行解码。
                                                                        */
/************************************************************************/
int DecodeIFrame(int16 encodeLumaData[XX][YY],FrameParament &curFrameData)
{
	int i,j;
	
	//将encodeLumaData数据写入sideLumaData中,作为下一帧解码的边信息
	memcpy(sideLumaData,encodeLumaData,XX*YY*sizeof(int16));

	for(i=0;i<X;i++)
	{
		for(j=0;j<Y;j++)
		{
			GetDecodeBlockData(i, j, encodeLumaData);
			SaveDecodeFile(i,j);
			DecodeMacroBlock();
			SetDecodeData(i,j,curFrameData); //将解码后的宏块保存到curFrameData.LumaData[XX][YY]
		}
	}
	return 0;
}

/************************************************************************/
/* 
   对P帧数据进行解码。
                                                                        */
/************************************************************************/
int DecodePFrame(int index,FrameParament &curFrameData,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS])
{
	int16 encodeLumaData[XX][YY] = {0};//用于存放当前帧未进行反量化反DCT操作的数据。
	int i,j;

	for(i=0;i<X;i++)
	{
		for(j=0;j<Y;j++)
		{
			GetSyndromeBits(i,j,syndromeLumaBits);
			GetBestSideInfo(i,j,syndromeCheck);//先通过Hash校验得出边信息再进行LDPC解码。
			theoryPb = 0.362296;// theoryPbFrame[index][i][j]; //0.183233;//
			SyndromeDecode3();
			SetQuantBlockData(i,j,uPrismLumaData,encodeLumaData);
			SaveDecodeFile(i,j);
			DecodeMacroBlock();
			SetDecodeData(i,j,curFrameData); //将解码后的宏块保存到curFrameData.LumaData[XX][YY]
		}
		fprintf(fp_mv,"\n");
	}

	//将当前帧解码出的LumaData数据写入sideLumaData中,作为下一帧解码的边信息
	memcpy(sideLumaData,encodeLumaData,XX*YY*sizeof(int16));

	return 0;
}

/************************************************************************/
/*
	获取一个用01序列表示的伴随子宏块序列。
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
	得到最佳伴随子信息
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
			if (SAD(x,y,i,j,syndromeCheck)==1) //通过SAD来选择边信息
			//if (SHD(x,y,i,j,syndromeCheck)==1)	 //通过汉明距来选择边信息
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
	
	//搜索的宏块区域超出一帧大小的边界。
	if( rx<0 || ry<0 || (rx+BLOCK_HEIGHT)>XX || (ry+BLOCK_WIDTH)>YY)
	{
//		fprintf(fp_sideInfo,"i=%d,j=%d,error!\n",x,y);
		return 0;
	}
	//得到以(rx,ry)为左顶点坐标的边信息宏块的十进制表示
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


	//将sideInfoSource数组里面的数转换成16bits的二进制数并存到数组sideInfo中
	tempIndex = 0;
	count = 0;
	for (i=0;i<PRISMBITS;i++)
	{
//		printf("%d\n",sideInfoSource[i]);
		tempIndex = count<<4;		//通过引入变量替代下面语句，减少计算。
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
	对一帧数据中的宏块进行伴随子解码,边信息在sideLumaData[XX][YY]数组中。
	伴随子解码后的数据存放在syndromeSource数组。
	                                                                    */
/************************************************************************/
int SyndromeDecode3()
{
	bool test =false;
	int count = 0;
	int tempIndex = 0;

	LDPCDecode();

	/*保存LDPC解码还原的二进制序列*/
	for (int k=0;k<cols;k++)
	{
		fprintf(fp_reSource,"%d ",decodeBlockBits[k]);
	}
	fprintf(fp_reSource,"\n");

	tempIndex = 0;
	count = 0;
	for (int i=0;i<PRISMBITS;i++)
	{
		tempIndex = count<<4;		//通过引入变量替代下面语句，减少计算。
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
	P帧中设置QuantBlockData和encodeLumaData数组
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
	LDPC解码部分。
	                                                                    */
/************************************************************************/
int LDPCDecode()
{
	int k=0;
	//通过伴随子序列syndromeBlockBits与边信息sideInfo进行解码。
	int iter;
	bool sign = false;
	uint8 xxSyndrome[rows] = {0};
	double q0[cols] = {0};
	double qin[cols][maxVariableDegree] = {0};
	double tin[rows][maxCheckDegree] = {0};
	
                   
	/*初始化q0*/
	double tempq0=log((1-theoryPb)/theoryPb)/log(2);
	for(k=0;k<cols;k++)
	{
		q0[k]=(1-2*sideInfo[k])*tempq0;
	}
	
	/*循环译码；*/
	for(iter=0;iter<MAX_ITER;iter++)
	{
//		printf("i=%d,j=%d,x=%d,y=%d,iter=%d,residual=%d,tempResidual=%d\n",i,j,x,y,iter,residual,tempRisual);

// 		fprintf(DecodeInfo,"i=%d,j=%d,x=%d,y=%d,iter=%d,residual=%d\n",i,j,x,y,iter,residual);
		qin2tin(tin,qin,q0);
		tin2qin(qin,tin,syndromeBlockBits);
	
		/*判决;*/
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
	
		if(sign)   /*sign用于标记x与xx是否相等,通过伴随子来判断;true表示完全正确译码;*/
		{
//			printf("iter=%d,正确译码\n",iter);
			finish = true;
			return 1;
		}
	}
	
	return 1;
}




/************************************************************************/
/*
	获取一个已经编码的宏块准备解码
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
	对一个宏块块进行解码，将解码反量化和反DCT结果存到MacroBlockData中
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
	对4x4子块进行反量化与反整数DCT变换
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

	qp_per    = (QP-MIN_QP)/6;	//被6除的整数部分
	qp_rem    = (QP-MIN_QP)%6;  //被6除的余数部分

	
	//反量化过程
	for (i=0;i<SUBBLOCK_HEIGHT;i++)
	{
		for (j=0;j<SUBBLOCK_WIDTH;j++)
		{
			W[i][j] = source[i][j]*dequant_coef[qp_rem][i][j]*(1<<qp_per);
		}
	}
	
	//将二维DCT反转换为两个一维变换先进行列变换
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
	
	//行变换
	for (i=0;i<SUBBLOCK_HEIGHT;i++)
    {
		temp[0][0] = xx[i][0] + xx[i][2];
		temp[0][1] = xx[i][0] - xx[i][2];
		temp[1][0] = xx[i][1] + (xx[i][3]>>1);
		temp[1][1] = (xx[i][1]>>1) - xx[i][3];

		//此处应该再除64，参考白皮书Transform&量化部分第五节解码部分。
		MacroBlockData[i+ox][0+oy] = (temp[0][0] + temp[1][0])>>6;   
		MacroBlockData[i+ox][1+oy] = (temp[0][1] + temp[1][1])>>6;
		MacroBlockData[i+ox][2+oy] = (temp[0][1] - temp[1][1])>>6;
		MacroBlockData[i+ox][3+oy] = (temp[0][0] - temp[1][0])>>6;
   }
	return 1;
}


/************************************************************************/
/* 
	解码后的宏块保存到curFrameData.LumaData[XX][YY]中;
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
	将一帧已经解码的数据存到DECODE_FILE文件中
                                                                        */
/************************************************************************/
int WriteOneDecodeFrame(int frameIndex,FILE *fp1,FrameParament curFrameData) 
{
	//存储ENCODE_FILE文件
	fseek(fp1,(XX*YY*3*frameIndex)>>1,SEEK_SET);
	if(fwrite(curFrameData.LumaData[0],sizeof(uint8),XX*YY,fp1)==0) return 0;
	if(fwrite(curFrameData.UData[0],sizeof(uint8),(XX*YY)>>2,fp1)==0) return 0;
	if(fwrite(curFrameData.VData[0],sizeof(uint8),(XX*YY)>>2,fp1)==0) return 0;
	return 1;
}


/************************************************************************/
/* 
	LDPC解码要用的函数
                                                                        */
/************************************************************************/
void qin2tin(double tin[rows][maxCheckDegree],double qin[cols][maxVariableDegree],double q0[cols])
{
	int i,j,k,l,m;
	double sum=0;
	double **qout;
/*	double qout[cols][maxVariableDegree]={0};*/
	bool flag=true;             /*设置标志位用于Operator是加还是乘;true对应"加"，false对应乘；*/
	double *temp;
	int *tempVariable;

	/************************************************************************/
	/* 为qout分配空间                                                       */
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

	for(i=0;i<cols;i++)    /*由qin得到qout*/
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
	for(k=0;k<rows;k++)		/*由qout得到tin	*/
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
				/*得到tin[i][j]中的i在qout中对应的序号；*/
				if(variableNodeOnes[variableNodeIndex][m]==k)
				{
					break;
				}
			}

			tin[k][l]=qout[variableNodeIndex][m];
		}

	}

	/************************************************************************/
	/* 释放动态分配的内存                                                   */
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
	/* 初始化tout                                                           */
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
	/* 释放动态分配的内存                                                   */
	for (i=0;i<rows;i++)
	{
		free(tout[i]);
	}
	free(tout);
	/************************************************************************/
}


/************************************************************************/
/* 根据flag来选择相应的运算                                             */
/************************************************************************/
double Operator(double *temp, int *temp1,int j,bool flag)
{
	int k;
	double tempSum;
	
	tempSum=0;

	if(flag)		/*运算符为加*/
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

/*运算符为乘;通过取对数函数变乘为加:ln(a*b)=ln(a)+ln(b)
			tempSum+=log(tanh((*(temp+k))/2));	
*/		}
		
	}
	return tempSum;
}


/***************************************************************************/
/* CompareSyndrome()用于比较两个伴随之是否相等，相等返回true，否则返回false*/
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
/* GetSyndrome()计算伴随子,并保存在syndrome指针所指的区域               */
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
		m = k>>4;		//k除BLOCK_WIDTH的整数部分
		n = k - (m<<4); //k除BLOCK_WIDTH的余数部分
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
		m = k>>4;		//k除BLOCK_WIDTH的整数部分
		n = k - (m<<4); //k除BLOCK_WIDTH的余数部分

		tempSad+=checkHamDistance(syndromeCheck[x][y*CHECKBITS+k],sideLumaData[rx+m][ry+n]);
	
	}
	if(tempSad >= sad)
		return 0;
	sad = tempSad;
	return 1;
}