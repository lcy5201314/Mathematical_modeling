 /* c7-4.h ����ͼ���ڽӶ��ر�洢��ʾ */
 #define MAX_VERTEX_NUM 20
 typedef enum{unvisited,visited}VisitIf;
 typedef struct EBox
 {
   VisitIf mark; /* ���ʱ�� */
   int ivex,jvex; /* �ñ����������������λ�� */
   struct EBox *ilink,*jlink; /* �ֱ�ָ�������������������һ���� */
   InfoType *info; /* �ñ���Ϣָ�� */
 }EBox;
 typedef struct
 {
   VertexType data;
   EBox *firstedge; /* ָ���һ�������ö���ı� */
 }VexBox;
 typedef struct
 {
   VexBox adjmulist[MAX_VERTEX_NUM];
   int vexnum,edgenum; /* ����ͼ�ĵ�ǰ�������ͱ��� */
 }AMLGraph;