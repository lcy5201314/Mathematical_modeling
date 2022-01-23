 /* c7-2.h ͼ���ڽӱ�洢��ʾ */
 #define MAX_VERTEX_NUM 20
 typedef enum{DG,DN,AG,AN}GraphKind; /* {����ͼ,������,����ͼ,������} */
 typedef struct ArcNode
 {
   int adjvex; /* �û���ָ��Ķ����λ�� */
   struct ArcNode *nextarc; /* ָ����һ������ָ�� */
   InfoType *info; /* ����Ȩֵָ�룩 */
 }ArcNode; /* ���� */
 typedef struct
 {
   VertexType data; /* ������Ϣ */
   ArcNode *firstarc; /* ��һ������ĵ�ַ,ָ���һ�������ö���Ļ���ָ�� */
 }VNode,AdjList[MAX_VERTEX_NUM]; /* ͷ��� */
 typedef struct
 {
   AdjList vertices;
   int vexnum,arcnum; /* ͼ�ĵ�ǰ�������ͻ��� */
   int kind; /* ͼ�������־ */
 }ALGraph;
