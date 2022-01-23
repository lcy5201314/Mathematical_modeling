 /* c7-3.h ����ͼ��ʮ������洢��ʾ */
 #define MAX_VERTEX_NUM 20
 typedef struct ArcBox
 {
   int tailvex,headvex; /* �û���β��ͷ�����λ�� */
   struct ArcBox *hlink,*tlink; /* �ֱ�Ϊ��ͷ��ͬ�ͻ�β��ͬ�Ļ������� */
   InfoType *info; /* �û������Ϣ��ָ��(����) */
 }ArcBox; /* ����� */
 typedef struct
 {
   VertexType data;
   ArcBox *firstin,*firstout; /* �ֱ�ָ��ö����һ���뻡�ͳ��� */
 }VexNode; /* ������ */
 typedef struct
 {
   VexNode xlist[MAX_VERTEX_NUM]; /* ��ͷ����(����) */
   int vexnum,arcnum; /* ����ͼ�ĵ�ǰ�������ͻ��� */
 }OLGraph;