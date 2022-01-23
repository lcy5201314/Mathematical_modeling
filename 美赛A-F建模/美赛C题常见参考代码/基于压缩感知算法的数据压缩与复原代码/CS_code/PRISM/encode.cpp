#include "encode.h"

uint8 MacroBlockData[BLOCK_HEIGHT][BLOCK_WIDTH];	//原始序列中的宏块
int16 QuantBlockData[BLOCK_HEIGHT][BLOCK_WIDTH];	//原始序列宏块经过DCT和量化
int16 DCTBlockData[BLOCK_HEIGHT][BLOCK_WIDTH];		//原始序列宏块经过DCT

int16 syndromeSource[PRISMBITS];					//要进行伴随子编码的宏块数据
int16 BlockSource[UPRISMBITS];						//要进行传统编码的宏块数据
uint8 syndromeBlockBits[rows];						//整个宏块的伴随字序列。
//uint8 syndromeLumaBits[X][Y][rows];				    //整个帧的伴随子序列。
//uint8 syndromeCheck[X][Y*CHECKBITS];			    //用于判断是否采用正确解码的Hash值。

/************************************************************************/
/*
	获取指定帧内的数据
					                                                    */
/************************************************************************/
int GetOneFrameData(int frameIndex,FILE *fp,FrameParament &curFrameData)
{
	fseek(fp,(XX*YY*3*frameIndex)>>1,SEEK_SET);
	if(fread(curFrameData.LumaData,sizeof(uint8),XX*YY,fp)==0) return 0;
	if(fread(curFrameData.UData,sizeof(uint8),(XX*YY)>>2,fp)==0) return 0;
	if(fread(curFrameData.VData,sizeof(uint8),(XX*YY)>>2,fp)==0) return 0;
	return 1;
}


/************************************************************************/
/* 
   对I帧数据进行编码。
                                                                        */
/************************************************************************/
int EncodeIFrame(FrameParament curFrameData,int16 encodeLumaData[XX][YY])
{
	int i,j;
	for(i=0;i<X;i++)
	{
		for(j=0;j<Y;j++)
		{
			GetBlockData(i,j,curFrameData);
			EncodeMacroBlock();
			SetEncodeData(i,j,encodeLumaData);			//将量化后的宏块保存到encodeLuma[XX*YY]
			SaveEncodeFile(i,j);
		}
	}
	return 1;
}

/************************************************************************/
/* 
   对P帧数据进行编码。
                                                                        */
/************************************************************************/
int EncodePFrame(FrameParament curFrameData,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS])
{
	int i,j;
	for(i=0;i<X;i++)
	{
		for(j=0;j<Y;j++)
		{
			GetBlockData(i,j,curFrameData);
			EncodeMacroBlock();
			GetSyndromeSource();
			SyndromeEncode();
			SetPrismData(i,j,uPrismLumaData,syndromeLumaBits,syndromeCheck);
			SaveEncodeFile(i,j);
		}
	}
	return 1;

}

/************************************************************************/
/*
	获取一个要进行编码的宏块数据
					                                                    */
/************************************************************************/
int GetBlockData(int x,int y,FrameParament curFrameData)
{	
	int ox = x<<4;                //i*BLOCK_HEIGHT;
	int oy = y<<4;			      //j*BLOCK_WIDTH; 
	for (int i=0;i<BLOCK_HEIGHT;i++)
	{
		for (int j=0;j<BLOCK_WIDTH;j++)
		{
			MacroBlockData[i][j] = curFrameData.LumaData[ox+i][oy+j];
		}
	}
	return 1;
}

/************************************************************************/
/*
	对I帧宏块进行编码
					                                                    */
/************************************************************************/
int EncodeMacroBlock()
{
	int i,j;
	uint8 subMacroBlock[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH] = {0};
	int16 subDctData[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH] = {0};
	for (i=0;i<subX;i++)
	{
		for (j=0;j<subY;j++)
		{
			for (int k=0;k<SUBBLOCK_HEIGHT;k++)
			{
				for (int l=0;l<SUBBLOCK_WIDTH;l++)
				{
					subMacroBlock[k][l] = MacroBlockData[(i<<2)+k][(j<<2)+l];
				}
			}
			DCT(subMacroBlock,subDctData);
			Quant(i,j,subDctData);
			//在此处加入zigzag扫描和熵编码;
		}
	}
	return 1;
}


/************************************************************************/
/* 
	对4x4子块进行整数DCT变换。
	                                                                    */
/************************************************************************/
int DCT(uint8 source[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH],int16 W[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH])
{
	int i,j,i2;
	int16 temp[2][2]={0};

	int16 xx[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH] = {0};

	//将二维DCT转换为两个一维变换先进行列变换
    for (i=0;i<SUBBLOCK_WIDTH;i++)
    {
		for (j=0;j<2;j++)
		{
			i2=3-j;
			temp[j][0] = source[j][i]+source[i2][i];
			temp[j][1] = source[j][i]-source[i2][i];
		}
		xx[0][i] = temp[0][0] + temp[1][0];
		xx[1][i] = (temp[0][1]<<1) + temp[1][1];
		xx[2][i] = temp[0][0] - temp[1][0];
		xx[3][i] = temp[0][1] -(temp[1][1]<<1);
    }
	
	//行变换
	for (i=0;i<SUBBLOCK_HEIGHT;i++)
    {
		for (j=0;j<2;j++)
		{
			i2=3-j;
			temp[j][0] = xx[i][j]+xx[i][i2];
			temp[j][1] = xx[i][j]-xx[i][i2];
		}
		W[i][0] = temp[0][0] + temp[1][0];
		W[i][1] = (temp[0][1]<<1) + temp[1][1];
		W[i][2] = temp[0][0] - temp[1][0];
		W[i][3] = temp[0][1] -(temp[1][1]<<1);
    }
	return 1;

}


/************************************************************************/
/* 
	对4x4子块进行量化；
	将DCT变换结果写入DCTBlockData[BLOCK_HEIGHT][BLOCK_WIDTH]；
	量化结果写入QuantBlockData[BLOCK_HEIGHT][BLOCK_WIDTH]。
	                                                                    */
/************************************************************************/
int Quant(int x,int y, int16 W[SUBBLOCK_HEIGHT][SUBBLOCK_WIDTH])
{
	int i,j;
	int ox = x<<2;         // x*SUBBLOCK_HEIGHT;
	int oy = y<<2;         // y*SUBBLOCK_WIDTH;
	int qp_per,qp_rem,q_bits,qp_const;
	qp_per = (QP-MIN_QP)/6;	//被6除的整数部分
	qp_rem = (QP-MIN_QP)%6;  //被6除的余数部分
	q_bits = Q_BITS+qp_per;
    qp_const=(1<<q_bits)/3;    // intra 对应白皮书中的f

	//量化过程
	for (i=0;i<SUBBLOCK_HEIGHT;i++)
	{
		for (j=0;j<SUBBLOCK_WIDTH;j++)
		{
			DCTBlockData[i+ox][j+oy] = W[i][j];
			QuantBlockData[i+ox][j+oy] = (abs(W[i][j])*quant_coef[qp_rem][i][j]+qp_const)>>q_bits;
			if (W[i][j]<0)
			{
				QuantBlockData[i+ox][j+oy] = -1*QuantBlockData[i+ox][j+oy];
			}
		}
	}
	return 1;
}



/************************************************************************/
/* 
	将I帧编码后宏块保存到encodeLuma[XX][YY];
	                                                                    */
/************************************************************************/
int SetEncodeData(int x,int y,int16 encodeLumaData[XX][YY])
{
	int ox,oy;
	int i,j;

	ox = x<<4;    //x*BLOCK_HEIGHT;
	oy = y<<4;    //y*BLOCK_WIDTH;
	for (i=0;i<BLOCK_HEIGHT;i++)
	{
		for (j=0;j<BLOCK_WIDTH;j++)
		{
			encodeLumaData[i+ox][j+oy] = QuantBlockData[i][j];
		}
	}

	return 1;
}

/************************************************************************/
/* 
	保存宏块进行DCT变换后的数据;
	                                                                    */
/************************************************************************/
void SaveEncodeFile(int x,int y)
{

	int i,j;
	fprintf(fp_enDCT,"i=%d,j=%d\n",x,y);
	fprintf(fp_enQuant,"i=%d,j=%d\n",x,y);


	for (i=0;i<BLOCK_HEIGHT;i++)
	{
		for (j=0;j<BLOCK_WIDTH;j++)
		{
			fprintf(fp_enDCT,"%4d ",DCTBlockData[i][j]);
			fprintf(fp_enQuant,"%4d ",QuantBlockData[i][j]);
		}
		fprintf(fp_enDCT,"\n");
		fprintf(fp_enQuant,"\n");
	}
}


/************************************************************************/
/* 
	分离要进行Wyner-ziv编码的数据和进行传统编码的数据
	                                                                    */
/************************************************************************/
int GetSyndromeSource()
{
	int tempIndex = 0;
	int synsourceIndex = 0;
	int blockSourceIndex = 0;
	int ox = 0;
	int oy = 0;

	tempIndex = 0;
	for (int i=0;i<subX;i++)
	{
		ox = i<<2;    //i*SUBBLOCK_HEIGHT;
		for (int j=0;j<subY;j++)
		{
			oy = j<<2;    //j*SUBBLOCK_WIDTH;
			synsourceIndex = tempIndex<<2;
			blockSourceIndex = tempIndex*12 - 4;
			syndromeSource[synsourceIndex+0] = QuantBlockData[ox+0][oy+0];
			syndromeSource[synsourceIndex+1] = QuantBlockData[ox+0][oy+1];
			syndromeSource[synsourceIndex+2] = QuantBlockData[ox+1][oy+0];
			syndromeSource[synsourceIndex+3] = QuantBlockData[ox+2][oy+0];
			BlockSource[blockSourceIndex+4] = QuantBlockData[ox+1][oy+1];
			BlockSource[blockSourceIndex+5] = QuantBlockData[ox+0][oy+2];
			BlockSource[blockSourceIndex+6] = QuantBlockData[ox+0][oy+3];
			BlockSource[blockSourceIndex+7] = QuantBlockData[ox+1][oy+2];
			BlockSource[blockSourceIndex+8] = QuantBlockData[ox+2][oy+1];
			BlockSource[blockSourceIndex+9] = QuantBlockData[ox+3][oy+0];
			BlockSource[blockSourceIndex+10] = QuantBlockData[ox+3][oy+1];
			BlockSource[blockSourceIndex+11] = QuantBlockData[ox+2][oy+2];
			BlockSource[blockSourceIndex+12] = QuantBlockData[ox+1][oy+3];
			BlockSource[blockSourceIndex+13] = QuantBlockData[ox+2][oy+3];
			BlockSource[blockSourceIndex+14] = QuantBlockData[ox+3][oy+2];
			BlockSource[blockSourceIndex+15] = QuantBlockData[ox+3][oy+3];
		
			tempIndex++;
		}
	}

	return 1;
}


/************************************************************************/
/* 
	设置hash校验子和伴随子等。
	                                                                    */
/************************************************************************/
int SetPrismData(int x,int y,int16 uPrismLumaData[X][Y][UPRISMBITS],
				 uint8 syndromeLumaBits[X][Y][rows],
				 int16 syndromeCheck[X][Y*CHECKBITS])
{
	for (int i=0;i<rows;i++)
	{
		syndromeLumaBits[x][y][i] = syndromeBlockBits[i];
	//	printf("%d ",syndromeBlockBits[i]);
	}
	//getch();
	for (int k=0;k<CHECKBITS;k++)
	{
		int m,n;
		m = k>>4;		//k除BLOCK_WIDTH的整数部分
		n = k - (m<<4); //k除BLOCK_WIDTH的余数部分
		//syndromeCheck[x][y*CHECKBITS+k] = DCTBlockData[m][n];	//设置CHECKBITS个校验位。
		syndromeCheck[x][y*CHECKBITS+k] = QuantBlockData[m][n];	//设置CHECKBITS个校验位。
	}

	for (i=0;i<UPRISMBITS;i++)
	{
		uPrismLumaData[x][y][i] = BlockSource[i];
	}

	return 1;
}

/************************************************************************/
/* 
	将宏块编成伴随子形式，宏块的伴随子序列保存在syndromeBlockBits
	                                                                    */
/************************************************************************/
int SyndromeEncode()
{
	int i,j;
	uint8 sourceBits[cols] = {0};  //每个像素点值换成16位的2进制

// 		for(i=0;i<rows;i++)	//检验是否得到校验矩阵
// 		{
// 			for(j=0;j<maxCheckDegree;j++)
// 			{
// 				printf("%d ",checkNodeOnes[i][j]);
// 			}
// 			getch();
// 		}

	int tempIndex = 0;
	int count =0;

	//将syndromeSource数组里面的数转换成16bits的二进制数

	for (i=0;i<PRISMBITS;i++)
	{
		tempIndex = (count++)<<4;		//通过引入变量替代下面语句，减少计算。
		//printf("%d\n",syndromeSource[i]);
		for (int k=0;k<16;k++)
		{
			if (syndromeSource[i] & 1<<(15-k))
			{
				sourceBits[tempIndex+k] = 1;
			}
			else
			{
				sourceBits[tempIndex+k] = 0;
			}
 		//		printf("%d ",sourceBits[tempIndex+k]);
 		//		getch();
		}
		//printf("\n");
	}

// 	sourceBits[cols-1] = 0; //补足1026位
// 	sourceBits[cols-2] = 0;

	/*保存要进行LDPC编码的原始二进制序列*/
	for (i=0;i<cols;i++)
	{
		fprintf(fp_source,"%d ",sourceBits[i]);
	}
	fprintf(fp_source,"\n");


	int tempSyndrome = 0;

	//计算伴随子序列
	for(i=0;i<rows;i++)	
	{
		tempSyndrome=0;
		for(j=0;j<maxCheckDegree;j++)
		{
			if(checkNodeOnes[i][j]==-1)
			{
				break;
			}
// 			printf("%d,%d ",checkNodeOnes[i][j],tempSyndrome);
// 			getch();
			tempSyndrome+=sourceBits[checkNodeOnes[i][j]];
		}
		syndromeBlockBits[i]=tempSyndrome%2;
	}

	return 1;
}

