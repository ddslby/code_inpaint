#pragma once


//���ظ�����
#ifndef INPAINTING_H
#define INPAINTING_H

#include "image.h"
#define save_path ".\\result"//���ر����޸�ͼ��

#define PAINT_COLOR  RGB(0,255,0)  //Ŀ��������ɫ(��ɫ)
#define SOURCE 0	//��õ�
#define TARGET -1	//�����
#define BOUNDARY -2	//�߽��

#define winsize 4	//���ڴ�С(4�ǰ뾶��ʵ�ʴ���Ϊ9*9)

typedef struct
{
	double grad_x;
	double grad_y;
}gradient; //�ݶȵĽṹ��

typedef struct
{
	double norm_x;
	double norm_y;
}norm;  // the structure that record the norm

//inpainting��
class inpainting
{
public://
	
	CImage * m_pImage;	//���޸�ͼ��
	bool m_bOpen;		//1��ͼ��ɹ�	0��ͼ��ʧ��
	int m_width;		//ͼ����
	int m_height;		//ͼ��߶�

	COLORREF * m_color;	//��ɫͼ��
	double * m_r;		//Rͨ��
	double * m_g;		//Gͨ��
	double * m_b;		//Bͨ��
	
	int m_top;			//�ϱ߽�
	int m_bottom;		//�±߽�
	int m_left;			//��߽�
	int m_right;		//�ұ߽�

	int * m_mark;		//SOURCE(0:��õ�)		TARGET(-1:�����)		BOUNDARY(-2:�߽��)
	double * m_confid;	//�������Ŷ�
	double * m_pri;		//�������ȼ�(���߽����ػ��õ�)
	double * m_gray;	//�Ҷ�ͼ��
	bool * m_source;	//�Ƿ������Ϊƥ���(0:������	1:����)

	inpainting(char * name);	//���캯��������ͼ������
	~inpainting(void);			//��������
	bool process(void);			//�����޸��㷨����Ҫʵ�ֹ���
	void DrawBoundary(void);	//ͼ���ϻ�����������߽�

	double ComputeConfidence(int i, int j);	//����(i,j)�����Ŷ�
	double priority(int x, int y);			//����(x,y)�����ȼ�
	double ComputeData(int i, int j);		//����(i,j)��������
	void Convert2Gray(void);				//������ͼ��ת��Ϊ�Ҷ�ͼ��
	gradient GetGradient(int i, int j);		//����(i,j)���ݶ�
	norm GetNorm(int i, int j);				// calculate the norm at one pixel��λ������
	bool draw_source(void);					//�ҳ�����ƥ���
	bool PatchTexture(int x, int y,int &patch_x,int &patch_y);//��������ҵ���(x,y)�������Ƶ��޲����λ��(patch_x,patch_y)
	bool update(int target_x, int target_y, int source_x, int source_y, double confid);	//�޸��˿鲢�����������Ŷ�
	bool TargetExist(void);					//�鿴�Ƿ�����Ҫ�޸���ͼ���
	void UpdateBoundary(int i, int j);		//���±߽���
	void UpdatePri(int i, int j);			//���±߽������ص����ȼ�
};


#endif
