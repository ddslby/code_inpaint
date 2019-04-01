#pragma once


//防重复定义
#ifndef INPAINTING_H
#define INPAINTING_H

#include "image.h"
#define save_path ".\\result"//本地保存修复图像

#define PAINT_COLOR  RGB(0,255,0)  //目标区域颜色(绿色)
#define SOURCE 0	//完好点
#define TARGET -1	//破损点
#define BOUNDARY -2	//边界点

#define winsize 4	//窗口大小(4是半径，实际窗口为9*9)

typedef struct
{
	double grad_x;
	double grad_y;
}gradient; //梯度的结构体

typedef struct
{
	double norm_x;
	double norm_y;
}norm;  // the structure that record the norm

//inpainting类
class inpainting
{
public://
	
	CImage * m_pImage;	//待修复图像
	bool m_bOpen;		//1打开图像成功	0打开图像失败
	int m_width;		//图像宽度
	int m_height;		//图像高度

	COLORREF * m_color;	//彩色图像
	double * m_r;		//R通道
	double * m_g;		//G通道
	double * m_b;		//B通道
	
	int m_top;			//上边界
	int m_bottom;		//下边界
	int m_left;			//左边界
	int m_right;		//右边界

	int * m_mark;		//SOURCE(0:完好点)		TARGET(-1:破损点)		BOUNDARY(-2:边界点)
	double * m_confid;	//像素置信度
	double * m_pri;		//像素优先级(仅边界像素会用到)
	double * m_gray;	//灰度图像
	bool * m_source;	//是否可以作为匹配块(0:不可以	1:可以)

	inpainting(char * name);	//构造函数，传入图像名称
	~inpainting(void);			//析构函数
	bool process(void);			//整个修复算法的主要实现过程
	void DrawBoundary(void);	//图像上绘制破损区域边界

	double ComputeConfidence(int i, int j);	//计算(i,j)的置信度
	double priority(int x, int y);			//计算(x,y)的优先级
	double ComputeData(int i, int j);		//计算(i,j)的数据项
	void Convert2Gray(void);				//将输入图像转换为灰度图像
	gradient GetGradient(int i, int j);		//计算(i,j)的梯度
	norm GetNorm(int i, int j);				// calculate the norm at one pixel单位法向量
	bool draw_source(void);					//找出所有匹配块
	bool PatchTexture(int x, int y,int &patch_x,int &patch_y);//从完好区找到与(x,y)块最相似的修补块的位置(patch_x,patch_y)
	bool update(int target_x, int target_y, int source_x, int source_y, double confid);	//修复此块并更新像素置信度
	bool TargetExist(void);					//查看是否还有需要修复的图像块
	void UpdateBoundary(int i, int j);		//更新边界线
	void UpdatePri(int i, int j);			//更新边界线像素的优先级
};


#endif
